#include "RandomDecimation.h"
#include "Utils.h"
#include "Test.h"

#include <iterator>
#include <limits>
#include <algorithm>
#include <chrono>

#include "pmp/algorithms/DistancePointTriangle.h"
#include "pmp/algorithms/SurfaceNormals.h"

using namespace pmp;
using namespace std::chrono;

namespace neuralSubdiv {

    RandomDecimation::RandomDecimation(SurfaceMesh& mesh)
        : mesh_(mesh), initialized_(false), queue_(nullptr), eng_(std::chrono::steady_clock::now().time_since_epoch().count())

    {
        aspect_ratio_ = 0;
        edge_length_ = 0;
        max_valence_ = 0;
        normal_deviation_ = 0;
        hausdorff_error_ = 0;

        // add properties
        vquadric_ = mesh_.add_vertex_property<Quadric>("v:quadric");
        eselected_ = mesh_.add_edge_property<bool>("e:selected");
        eskip_ = mesh_.add_halfedge_property<bool>("he:skip");

        // get properties
        vpoint_ = mesh_.vertex_property<Point>("v:point");

        // compute face normals
        SurfaceNormals::compute_face_normals(mesh_);
        fnormal_ = mesh_.face_property<Normal>("f:normal");
    }

    RandomDecimation::~RandomDecimation()
    {
        // remove added properties
        mesh_.remove_vertex_property(vquadric_);
        mesh_.remove_face_property(normal_cone_);
        mesh_.remove_face_property(face_points_);
        mesh_.remove_edge_property(eselected_);
        mesh_.remove_halfedge_property(eskip_);
    }

    void RandomDecimation::initialize(Scalar aspect_ratio, Scalar edge_length,
        unsigned int max_valence,
        Scalar normal_deviation,
        Scalar hausdorff_error,
        unsigned int edge_subset_size,
        bool use_subset)
    {
        if (!mesh_.is_triangle_mesh())
            return;

        // store parameters
        aspect_ratio_ = aspect_ratio;
        max_valence_ = max_valence;
        edge_length_ = edge_length;
        normal_deviation_ = normal_deviation / 180.0 * M_PI;
        hausdorff_error_ = hausdorff_error;
        edge_subset_size_ = edge_subset_size;
        use_subset_ = use_subset;

        // properties
        if (normal_deviation_ > 0.0)
            normal_cone_ = mesh_.face_property<NormalCone>("f:normalCone");
        else
            mesh_.remove_face_property(normal_cone_);
        if (hausdorff_error > 0.0)
            face_points_ = mesh_.face_property<Points>("f:points");
        else
            mesh_.remove_face_property(face_points_);

        // vertex selection
        vselected_ = mesh_.get_vertex_property<bool>("v:selected");
        has_selection_ = false;
        if (vselected_)
        {
            for (auto v : mesh_.vertices())
            {
                if (vselected_[v])
                {
                    has_selection_ = true;
                    break;
                }
            }
        }

        // feature vertices/edges
        has_features_ = false;
        vfeature_ = mesh_.get_vertex_property<bool>("v:feature");
        efeature_ = mesh_.get_edge_property<bool>("e:feature");
        if (vfeature_ && efeature_)
        {
            for (auto v : mesh_.vertices())
            {
                if (vfeature_[v])
                {
                    has_features_ = true;
                    break;
                }
            }
        }

        // initialize quadrics
        for (auto v : mesh_.vertices())
        {
            vquadric_[v].clear();

            if (!mesh_.is_isolated(v))
            {
                for (auto f : mesh_.faces(v))
                {
                    vquadric_[v] += Quadric(fnormal_[f], vpoint_[v]);
                }
            }
        }

        // initialize normal cones
        if (normal_deviation_)
        {
            for (auto f : mesh_.faces())
            {
                normal_cone_[f] = NormalCone(fnormal_[f]);
            }
        }

        // initialize faces' point list
        if (hausdorff_error_)
        {
            for (auto f : mesh_.faces())
            {
                Points().swap(face_points_[f]); // free mem
            }
        }

        initialized_ = true;
    }

    void RandomDecimation::simplify(unsigned int n_vertices)
    {
        if (n_vertices >= mesh_.n_vertices())
        {
            std::cerr << "Simplification error: target number of vertices must be inferior to the current number of vertices." << std::endl;
            return;
        }

        if (!mesh_.is_triangle_mesh())
        {
            std::cerr << "Not a triangle mesh!" << std::endl;
            return;
        }

        // make sure the decimater is initialized
        if (!initialized_)
            initialize();

        unsigned int nv(mesh_.n_vertices());

        std::vector<Vertex> one_ring;
        std::vector<Vertex>::iterator or_it, or_end;
        Halfedge h;
        Vertex v;
        int it = 0;

        int dec_info_size = mesh_.n_vertices() - n_vertices;
        dec_infos_.resize(dec_info_size);
#ifdef DEBUG_PRINT
        std::cout << "estimation: " << dec_info_size << std::endl;
#endif

        // add properties for priority queue
        vpriority_ = mesh_.add_vertex_property<float>("v:prio");
        heap_pos_ = mesh_.add_vertex_property<int>("v:heap");
        vtarget_ = mesh_.add_vertex_property<Halfedge>("v:target");

        // select random edge subset
        if (use_subset_)
        {
            if (mesh_.n_edges() < edge_subset_size_)
            {
                use_subset_ = false;
            }
            else 
            {
                int idx = 0;
                edges_copy.resize(mesh_.n_edges());
                for (auto eit : mesh_.edges())
                {
                    edges_copy[idx] = eit.idx();
                    ++idx;
                }
                std::shuffle(edges_copy.begin(), edges_copy.end(), eng_);
            }
        }

        // build priority queue
        Halfedge he1, he2;
        float min_prio = std::numeric_limits<float>::max();
        Halfedge min_h;
        pmp::SurfaceMesh save;
        int stall_count = 0;
        while (nv > n_vertices)
        {
            if (it%50==0)std::cout << "it=" << it << std::endl;
            auto start_loop = high_resolution_clock::now();
            stall_count++;
            // get 1st element
            /*v = queue_->front();
            queue_->pop_front();*/
            auto start_prio= high_resolution_clock::now();

            min_prio = std::numeric_limits<float>::max();
            min_h = pmp::Halfedge(-1);
            for (int i = 0; i < 100; ++i)
            {
                he1 = mesh_.halfedge(pmp::Edge(edges_copy[i]), 0);
                he2 = mesh_.halfedge(pmp::Edge(edges_copy[i]), 1);
                CollapseData cd1(mesh_, he1);
                CollapseData cd2(mesh_, he2);
                if (is_collapse_legal(cd1) && priority(cd1) < min_prio)
                {
                    min_prio = priority(cd1);
                    min_h = he1;
                }
                if (is_collapse_legal(cd2) && priority(cd2) < min_prio)
                {
                    min_prio = priority(cd2);
                    min_h = he2;
                }
            }
            if (min_h.idx() == -1)
            {
                std::cout << "WARNING: premature stop in decimation (no more decent collapse available)" << std::endl;
                break;
            }
            auto stop_prio = high_resolution_clock::now();

            CollapseData cd(mesh_, min_h);
            /*h = vtarget_[v];
            (mesh_, h);*/

            // check this (again)
            if (!mesh_.is_collapse_ok(min_h))
                continue;
            // store one-ring
            one_ring.clear();
            for (auto vv : mesh_.vertices(cd.v0))
                one_ring.push_back(vv);

#ifdef DEBUG_PRINT
            std::cout << nv << std::endl;
#endif

            // naming convention: "after" = post edge collapse
            // meaning "onering" now means the one ring of vertex vi 
            // instead of the one ring of the edge vi,vj
            //auto start_param = high_resolution_clock::now();

            Eigen::MatrixXi F(mesh_.n_faces(), 3);
            neuralSubdiv::faces_to_matrix(mesh_, F);

            pmp::Vertex vi = cd.v1; pmp::Vertex vj = cd.v0;     // vi=remaining vertex
            Eigen::Vector2i boundary_idx;                       // idx of vi, vj in uv
            Eigen::MatrixXi F_uv, F_uv_after;                   // face list (idx in uv)
            Eigen::MatrixXi F_onering, V_map, F_map;            // mapping idx uv -> idx in V
            Eigen::MatrixXd uv, uv_after;                       // one ring vertices position in 2D
            Eigen::MatrixXd boundary_constraints;               // 2D position of the constraints

            neuralSubdiv::flatten_one_ring(mesh_, vi, vj,
                uv, F_uv, F_onering,
                boundary_idx, boundary_constraints,
                V_map, F_map);

            //assert(neuralSubdiv::check_F_uv(mesh_, F_map, F_uv, V_map));

            //std::cout << "F" << std::endl;
            //std::cout << F << std::endl;

            //for (int j = 0; j < F_uv.rows(); ++j)
            //{
            //    std::cout << pmp::Vertex(F_uv(j, 0)) << std::endl;
            //    std::cout << pmp::Vertex(F_uv(j, 1)) << std::endl;
            //    std::cout << pmp::Vertex(F_uv(j, 2)) << std::endl;
            //    std::cout << mesh_.find_halfedge(pmp::Vertex(F_uv(j, 2)), pmp::Vertex(F_uv(j, 0))) << std::endl;
            //    std::cout << mesh_.find_halfedge(pmp::Vertex(F_uv(j, 1)), pmp::Vertex(F_uv(j, 0))) << std::endl;
            //    std::cout << mesh_.find_halfedge(pmp::Vertex(F_uv(j, 2)), pmp::Vertex(F_uv(j, 1))) << std::endl;
            //    /*std::cout << (mesh_.find_edge(pmp::Vertex(F_uv(j, 0)), pmp::Vertex(F_uv(j, 1)))) << std::endl;
            //    std::cout << (mesh_.find_edge(pmp::Vertex(F_uv(j, 1)), pmp::Vertex(F_uv(j, 0)))) << std::endl;
            //    std::cout << (mesh_.find_edge(pmp::Vertex(F_uv(j, 1)), pmp::Vertex(F_uv(j, 2)))) << std::endl;
            //    std::cout << (mesh_.find_edge(pmp::Vertex(F_uv(j, 2)), pmp::Vertex(F_uv(j, 1)))) << std::endl;
            //    std::cout << (mesh_.find_edge(pmp::Vertex(F_uv(j, 2)), pmp::Vertex(F_uv(j, 0)))) << std::endl;
            //    std::cout << (mesh_.find_edge(pmp::Vertex(F_uv(j, 0)), pmp::Vertex(F_uv(j, 2)))) << std::endl;*/
            //    assert((mesh_.find_edge(pmp::Vertex(F_uv(j, 0)), pmp::Vertex(F_uv(j, 1)))).is_valid()
            //    || (mesh_.find_edge(pmp::Vertex(F_uv(j, 1)), pmp::Vertex(F_uv(j, 0)))).is_valid());
            //    assert((mesh_.find_edge(pmp::Vertex(F_uv(j, 2)), pmp::Vertex(F_uv(j, 1)))).is_valid()
            //    || (mesh_.find_edge(pmp::Vertex(F_uv(j, 1)), pmp::Vertex(F_uv(j, 2)))).is_valid());
            //    assert((mesh_.find_edge(pmp::Vertex(F_uv(j, 2)), pmp::Vertex(F_uv(j, 0)))).is_valid(),
            //        || (mesh_.find_edge(pmp::Vertex(F_uv(j, 0)), pmp::Vertex(F_uv(j, 2)))).is_valid());
            //}
            
#ifdef DEBUG_PRINT
            std::cout << "uv" << std::endl;
            std::cout << uv << std::endl;
#endif 

            if (!neuralSubdiv::check_lscm_self_folding(uv, F_uv, boundary_idx))
            {
                std::cerr << "Self-folding check failed (1st lscm)" << std::endl;
            }

            save = pmp::SurfaceMesh(mesh_);
            save.collapse(min_h);

            Eigen::MatrixXi V_map_after(save.valence(vi) + 1, 1);
            Eigen::MatrixXi F_onering_after;
            flatten_one_ring_after(save, vi, uv, V_map, uv_after, F_uv_after, V_map_after);

            bool tri_quality = neuralSubdiv::check_triangle_quality(save, vi);
            if (!tri_quality)
            {
                std::cout << "TRI bad tri quality" << std::endl;
                eskip_[min_h] = true;
                if (stall_count > 500)
                {
                    std::cout << "WARNING: premature stop in decimation (no more decent collapse available)" << std::endl;
                    break;
                }
                continue;
            }

            bool uv_quality = neuralSubdiv::check_triangle_quality_uv(uv_after, F_uv_after);
            if (!uv_quality)
            {
                std::cout << "UV bad tri quality" << std::endl;
                eskip_[min_h] = true;
                if (stall_count > 500)
                {
                    std::cout << "WARNING: premature stop in decimation (no more decent collapse available)" << std::endl;
                    break;
                }
                continue;
            }

            mesh_.collapse(min_h);
            --nv;
            postprocess_collapse(cd); // postprocessing, e.g., update quadrics
            //std::cout << "---" << std::endl;
            //flatten_one_ring_after(mesh_, vi, uv, V_map, uv_after, F_uv_after, V_map_after);
            //tri_quality = neuralSubdiv::check_triangle_quality(mesh_, vi);
            //if (!tri_quality)
            //{
            //    std::cout << "TRI WTF" << std::endl;
            //    assert(false);
            //}

            //uv_quality = neuralSubdiv::check_triangle_quality_uv(uv_after, F_uv_after);
            //if (!uv_quality)
            //{
            //    std::cout << "UV WTF" << std::endl;
            //    assert(false);
            //}

            //assert(neuralSubdiv::check_F_uv_after(mesh_, cd.v1, F_uv_after, V_map_after));


            //auto stop_param = high_resolution_clock::now();

#ifdef DEBUG_PRINT
            std::cout << "uv after" << std::endl;
            std::cout << uv_after << std::endl;
#endif
            DecInfo info;
            info.n_collapse = it;
            info.vi = cd.v1;
            info.boundary_idx = boundary_idx;
            info.uv = uv;
            info.uv_after = uv_after;
            info.F_uv = F_uv;
            info.F_uv_after = F_uv_after;
            info.V_map = V_map;
            info.V_map_after = V_map_after;
            dec_infos_[it] = info;

            // clean previous edge selection and recompute selection
            //auto start_r = high_resolution_clock::now();

            // TODO at the beginning of the loop
            if (use_subset_)
            {
                //for (int i = 0; i < 100; ++i)
                //    eselected_[pmp::Edge(edges_copy[i])] = false;
                if (mesh_.n_edges() < edge_subset_size_)
                {
                    use_subset_ = false;
                }
                else
                {
                    int idx = 0;
                    edges_copy.resize(mesh_.n_edges());
                    for (auto eit : mesh_.edges())
                    {
                        edges_copy[idx] = eit.idx();
                        ++idx;
                    }
                    std::shuffle(edges_copy.begin(), edges_copy.end(), eng_);
                    //for (int i = 0; i < 100; ++i)
                    //    eselected_[pmp::Edge(edges_copy[i])] = true;
                    //queue_->clear();
                    //queue_->reserve(100);
                    //for (int i = 0; i < 100; ++i)
                    //{
                    //    queue_->reset_heap_position(mesh_.vertex(pmp::Edge(edges_copy[i]), 0));
                    //    enqueue_vertex(mesh_.vertex(pmp::Edge(edges_copy[i]), 0));
                    //    queue_->reset_heap_position(mesh_.vertex(pmp::Edge(edges_copy[i]), 1));
                    //    enqueue_vertex(mesh_.vertex(pmp::Edge(edges_copy[i]), 1));
                    //}
                    //std::cout << "queue_ size:" << queue_->size() << std::endl;
                }
            }
            //auto stop_r = high_resolution_clock::now();

            // update queue
            //for (or_it = one_ring.begin(), or_end = one_ring.end(); or_it != or_end;
            //    ++or_it)
            //{
            //    eselected_[mesh_.edge(mesh_.halfedge(*or_it))] = true;
            //    enqueue_vertex(*or_it);
            //}
            stall_count = 0;
            ++it; // n_collapse 
            //auto stop_loop = high_resolution_clock::now();
            //auto duration_loop = duration_cast<milliseconds>(stop_loop - start_loop);
            //auto duration_r = duration_cast<milliseconds>(stop_r - start_r);
            //auto duration_prio = duration_cast<microseconds>(stop_prio - start_prio);
            //auto duration_param = duration_cast<milliseconds>(stop_param - start_param);
            //auto duration_f= duration_cast<microseconds>(stop_f- start_f);
            //std::cout << "Loop: " << duration_loop.count() << std::endl;
            //std::cout << "rand: " << duration_r.count() << std::endl;
            //std::cout << "prio: " << duration_prio.count() << std::endl;
            //std::cout << "param: " << duration_param.count() << std::endl;
            //std::cout << "face: " << duration_f.count() << std::endl;
        }

        // clean up
        //delete queue_;

#ifdef DEBUG_PRINT
        std::cout << "dec info size : " << dec_infos.size() << std::endl;
#endif
    }

    void RandomDecimation::enqueue_vertex(Vertex v)
    {
        float prio, min_prio(std::numeric_limits<float>::max());
        Halfedge min_h;

        // find best out-going halfedge
        for (auto h : mesh_.halfedges(v))
        {
            CollapseData cd(mesh_, h);
            if (is_collapse_legal(cd))
            {
                prio = priority(cd);
                if (prio != -1.0 && prio < min_prio)
                {
                    min_prio = prio;
                    min_h = h;
                }
            }
        }

        // target found -> put vertex on heap
        if (min_h.is_valid())
        {
            vpriority_[v] = min_prio;
            vtarget_[v] = min_h;

            if (queue_->is_stored(v))
                queue_->update(v);
            else
                queue_->insert(v);
        }

        // not valid -> remove from heap
        else
        {
            if (queue_->is_stored(v))
                queue_->remove(v);

            vpriority_[v] = -1;
            vtarget_[v] = min_h;
        }
    }

    bool RandomDecimation::is_collapse_legal(const CollapseData& cd)
    {
        // test selected vertices
        if (has_selection_)
        {
            if (!vselected_[cd.v0])
                return false;
        }

        if (eskip_[cd.v0v1])
            return false;

        // test features
        if (has_features_)
        {
            if (vfeature_[cd.v0] && !efeature_[mesh_.edge(cd.v0v1)])
                return false;

            if (cd.vl.is_valid() && efeature_[mesh_.edge(cd.vlv0)])
                return false;

            if (cd.vr.is_valid() && efeature_[mesh_.edge(cd.v0vr)])
                return false;
        }

        // do not collapse boundary vertices to interior vertices
        if (mesh_.is_boundary(cd.v0) && !mesh_.is_boundary(cd.v1))
            return false;

        // there should be at least 2 incident faces at v0
        if (mesh_.cw_rotated_halfedge(mesh_.cw_rotated_halfedge(cd.v0v1)) ==
            cd.v0v1)
            return false;

        // topological check
        if (!mesh_.is_collapse_ok(cd.v0v1))
            return false;

        // check maximal valence
        if (max_valence_ > 0)
        {
            unsigned int val0 = mesh_.valence(cd.v0);
            unsigned int val1 = mesh_.valence(cd.v1);
            unsigned int val = val0 + val1 - 1;
            if (cd.fl.is_valid())
                --val;
            if (cd.fr.is_valid())
                --val;
            if (val > max_valence_ && !(val < std::max(val0, val1)))
                return false;
        }

        // remember the positions of the endpoints
        const Point p0 = vpoint_[cd.v0];
        const Point p1 = vpoint_[cd.v1];

        // check for maximum edge length
        if (edge_length_)
        {
            for (auto v : mesh_.vertices(cd.v0))
            {
                if (v != cd.v1 && v != cd.vl && v != cd.vr)
                {
                    if (norm(vpoint_[v] - p1) > edge_length_)
                        return false;
                }
            }
        }

        // check for flipping normals
        if (normal_deviation_ == 0.0)
        {
            vpoint_[cd.v0] = p1;
            for (auto f : mesh_.faces(cd.v0))
            {
                if (f != cd.fl && f != cd.fr)
                {
                    Normal n0 = fnormal_[f];
                    Normal n1 = SurfaceNormals::compute_face_normal(mesh_, f);
                    if (dot(n0, n1) < 0.0)
                    {
                        vpoint_[cd.v0] = p0;
                        return false;
                    }
                }
            }
            vpoint_[cd.v0] = p0;
        }

        // check normal cone
        else
        {
            vpoint_[cd.v0] = p1;

            Face fll, frr;
            if (cd.vl.is_valid())
                fll = mesh_.face(
                    mesh_.opposite_halfedge(mesh_.prev_halfedge(cd.v0v1)));
            if (cd.vr.is_valid())
                frr = mesh_.face(
                    mesh_.opposite_halfedge(mesh_.next_halfedge(cd.v1v0)));

            for (auto f : mesh_.faces(cd.v0))
            {
                if (f != cd.fl && f != cd.fr)
                {
                    NormalCone nc = normal_cone_[f];
                    nc.merge(SurfaceNormals::compute_face_normal(mesh_, f));

                    if (f == fll)
                        nc.merge(normal_cone_[cd.fl]);
                    if (f == frr)
                        nc.merge(normal_cone_[cd.fr]);

                    if (nc.angle() > 0.5 * normal_deviation_)
                    {
                        vpoint_[cd.v0] = p0;
                        return false;
                    }
                }
            }

            vpoint_[cd.v0] = p0;
        }

        // check aspect ratio
        if (aspect_ratio_)
        {
            Scalar ar0(0), ar1(0);

            for (auto f : mesh_.faces(cd.v0))
            {
                if (f != cd.fl && f != cd.fr)
                {
                    // worst aspect ratio after collapse
                    vpoint_[cd.v0] = p1;
                    ar1 = std::max(ar1, aspect_ratio(f));
                    // worst aspect ratio before collapse
                    vpoint_[cd.v0] = p0;
                    ar0 = std::max(ar0, aspect_ratio(f));
                }
            }

            // aspect ratio is too bad, and it does also not improve
            if (ar1 > aspect_ratio_ && ar1 > ar0)
                return false;
        }

        // check Hausdorff error
        if (hausdorff_error_)
        {
            Points points;
            bool ok;

            // collect points to be tested
            for (auto f : mesh_.faces(cd.v0))
            {
                std::copy(face_points_[f].begin(), face_points_[f].end(),
                    std::back_inserter(points));
            }
            points.push_back(vpoint_[cd.v0]);

            // test points against all faces
            vpoint_[cd.v0] = p1;
            for (auto point : points)
            {
                ok = false;

                for (auto f : mesh_.faces(cd.v0))
                {
                    if (f != cd.fl && f != cd.fr)
                    {
                        if (distance(f, point) < hausdorff_error_)
                        {
                            ok = true;
                            break;
                        }
                    }
                }

                if (!ok)
                {
                    vpoint_[cd.v0] = p0;
                    return false;
                }
            }
            vpoint_[cd.v0] = p0;
        }

        // collapse passed all tests -> ok
        return true;
    }

    float RandomDecimation::priority(const CollapseData& cd)
    {
        // computer quadric error metric
        Quadric Q = vquadric_[cd.v0];
        Q += vquadric_[cd.v1];
        return Q(vpoint_[cd.v1]);
    }

    void RandomDecimation::postprocess_collapse(const CollapseData& cd)
    {
        // update error quadrics
        vquadric_[cd.v1] += vquadric_[cd.v0];

        // update normal cones
        if (normal_deviation_)
        {
            for (auto f : mesh_.faces(cd.v1))
            {
                normal_cone_[f].merge(
                    SurfaceNormals::compute_face_normal(mesh_, f));
            }

            if (cd.vl.is_valid())
            {
                Face f = mesh_.face(cd.v1vl);
                if (f.is_valid())
                    normal_cone_[f].merge(normal_cone_[cd.fl]);
            }

            if (cd.vr.is_valid())
            {
                Face f = mesh_.face(cd.vrv1);
                if (f.is_valid())
                    normal_cone_[f].merge(normal_cone_[cd.fr]);
            }
        }

        // update Hausdorff error
        if (hausdorff_error_)
        {
            Points points;

            // collect points to be distributed

            // points of v1's one-ring
            for (auto f : mesh_.faces(cd.v1))
            {
                std::copy(face_points_[f].begin(), face_points_[f].end(),
                    std::back_inserter(points));
                face_points_[f].clear();
            }

            // points of the 2 removed triangles
            if (cd.fl.is_valid())
            {
                std::copy(face_points_[cd.fl].begin(), face_points_[cd.fl].end(),
                    std::back_inserter(points));
                Points().swap(face_points_[cd.fl]); // free mem
            }
            if (cd.fr.is_valid())
            {
                std::copy(face_points_[cd.fr].begin(), face_points_[cd.fr].end(),
                    std::back_inserter(points));
                Points().swap(face_points_[cd.fr]); // free mem
            }

            // the removed vertex
            points.push_back(vpoint_[cd.v0]);

            // test points against all faces
            Scalar d, dd;
            Face ff;

            for (auto point : points)
            {
                dd = std::numeric_limits<Scalar>::max();

                for (auto f : mesh_.faces(cd.v1))
                {
                    d = distance(f, point);
                    if (d < dd)
                    {
                        ff = f;
                        dd = d;
                    }
                }

                face_points_[ff].push_back(point);
            }
        }
    }

    Scalar RandomDecimation::aspect_ratio(Face f) const
    {
        // min height is area/maxLength
        // aspect ratio = length / height
        //              = length * length / area

        SurfaceMesh::VertexAroundFaceCirculator fvit = mesh_.vertices(f);

        const Point p0 = vpoint_[*fvit];
        const Point p1 = vpoint_[*(++fvit)];
        const Point p2 = vpoint_[*(++fvit)];

        const Point d0 = p0 - p1;
        const Point d1 = p1 - p2;
        const Point d2 = p2 - p0;

        const Scalar l0 = sqrnorm(d0);
        const Scalar l1 = sqrnorm(d1);
        const Scalar l2 = sqrnorm(d2);

        // max squared edge length
        const Scalar l = std::max(l0, std::max(l1, l2));

        // triangle area
        Scalar a = norm(cross(d0, d1));

        return l / a;
    }

    Scalar RandomDecimation::distance(Face f, const Point& p) const
    {
        SurfaceMesh::VertexAroundFaceCirculator fvit = mesh_.vertices(f);

        const Point p0 = vpoint_[*fvit];
        const Point p1 = vpoint_[*(++fvit)];
        const Point p2 = vpoint_[*(++fvit)];

        Point n;

        return dist_point_triangle(p, p0, p1, p2, n);
    }

    void RandomDecimation::get_random_edge_subset(int n, std::vector<int>& edges_subset) 
    {
        edges_subset.resize(mesh_.n_edges());
        int idx = 0;
        for (auto eit : mesh_.edges())
        {
            edges_subset[idx] = eit.idx();
            idx++;
        }
        std::random_shuffle(edges_copy.begin(), edges_copy.end());
    }

    RandomDecimation::CollapseData::CollapseData(SurfaceMesh& sm, Halfedge h)
        : mesh(sm)
    {
        v0v1 = h;
        v1v0 = mesh.opposite_halfedge(v0v1);
        v0 = mesh.to_vertex(v1v0);
        v1 = mesh.to_vertex(v0v1);
        fl = mesh.face(v0v1);
        fr = mesh.face(v1v0);

        // get vl
        if (fl.is_valid())
        {
            v1vl = mesh.next_halfedge(v0v1);
            vlv0 = mesh.next_halfedge(v1vl);
            vl = mesh.to_vertex(v1vl);
        }

        // get vr
        if (fr.is_valid())
        {
            v0vr = mesh.next_halfedge(v1v0);
            vrv1 = mesh.prev_halfedge(v0vr);
            vr = mesh.from_vertex(vrv1);
        }
    }

} // namespace neuralSubdiv
