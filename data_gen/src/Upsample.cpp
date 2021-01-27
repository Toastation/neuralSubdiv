#include "Upsample.h"
#include "Utils.h"

#include <igl/find.h>

namespace neuralSubdiv {



    // based on pmp::SurfaceSubdivision::loop
    void upsample_mid_point(pmp::SurfaceMesh& mesh, pmp::SurfaceMesh& coarse)
    {
        pmp::VertexProperty<pmp::Point>points_ = mesh.vertex_property<pmp::Point>("v:point");
        pmp::VertexProperty<bool> vfeature_ = mesh.get_vertex_property<bool>("v:feature");
        pmp::EdgeProperty<bool> efeature_ = mesh.get_edge_property<bool>("e:feature");
        pmp::VertexProperty<Eigen::Vector3i> bary_face_prop = mesh.get_vertex_property<Eigen::Vector3i>("v:bary_face");
        pmp::VertexProperty<Eigen::Vector3d> bary_coord_prop = mesh.get_vertex_property<Eigen::Vector3d>("v:bary_coord");
        pmp::VertexProperty<bool> oldvert_prop = mesh.get_vertex_property<bool>("v:oldvert");

        if (!mesh.is_triangle_mesh())
            return;

        // reserve memory
        size_t nv = mesh.n_vertices();
        size_t ne = mesh.n_edges();
        size_t nf = mesh.n_faces();
        mesh.reserve(nv + ne, 2 * ne + 3 * nf, 4 * nf);

        // add properties
        auto vpoint = mesh.add_vertex_property<pmp::Point>("loop:vpoint");
        auto epoint = mesh.add_edge_property<pmp::Point>("loop:epoint");
        bool interior_edge = false;

        for (auto e : mesh.edges())
            epoint[e] = (points_[mesh.vertex(e, 0)] + points_[mesh.vertex(e, 1)]) * pmp::Scalar(0.5);

        pmp::Halfedge e_he;
        pmp::Vertex v0(-1), v1(-1);
        bool old0, old1;
        Eigen::Vector3i bf0, bf1, bf_final;
        Eigen::Vector3d bco0, bco1, bco_final;

        // insert new vertices on edges
        for (auto e : mesh.edges())
        {
            e_he = mesh.halfedge(e, 0);
            v0 = mesh.from_vertex(e_he);
            v1 = mesh.to_vertex(e_he);
            old0 = oldvert_prop[v0]; old1 = oldvert_prop[v1];
            bf0 = bary_face_prop[v0]; bf1 = bary_face_prop[v1];
            bco0 = bary_coord_prop[v0]; bco1 = bary_coord_prop[v1];

            // compute new vertex bary coord on the coarse mesh
            // TODO: remove this if, i think this case is already handled below
            if (old0 && old1)
            {
                bf_final(0) = v0.idx();
                bf_final(1) = v1.idx();
                for (auto v_it : mesh.vertices(mesh.face(e_he)))
                    if (v_it != v0 && v_it != v1 && oldvert_prop[v_it])
                        bf_final(2) = v_it.idx();
                bco_final(0) = 0.5;
                bco_final(1) = 0.5;
                bco_final(2) = 0;
            }
            else
            {
                bool in_list = false;
                int idx = 0;

                // add non zero bary coord of v0
                for (int i = 0; i < 3; ++i)
                {
                    if (bco0(i) >= 1e-7)
                    {
                        bf_final(idx) = bf0(i);
                        bco_final(idx) = bco0(i);
                        ++idx;
                    }
                }

                // add non zero bary coord of v1
                for (int i = 0; i < 3; ++i)
                {
                    if (bco1(i) >= 1e-7)
                    {
                        // already exists? sum
                        in_list = false;
                        for (int j = 0; j < idx; ++j)
                        {
                            if (bf1(i) == bf_final(j))
                            {
                                bco_final(j) += bco1(i);
                                in_list = true;
                                break;
                            }
                        }
                        if (!in_list)
                        {
                            assert(idx <= 2);
                            bf_final(idx) = bf1(i);
                            bco_final(idx) = bco1(i);
                            ++idx;
                        }
                    }
                }
                assert(idx >= 2);
                
                // if there is no third vertex, find the one in an adjacent face in the coarse mesh
                if (idx == 2)
                {
                    auto he_tmp = coarse.find_halfedge(pmp::Vertex(bf_final(0)), pmp::Vertex(bf_final(1)));
                    assert(he_tmp.is_valid());
                    for (auto v_it : coarse.vertices(coarse.face(he_tmp)))
                    {
                        if (v_it.idx() != bf_final(0) && v_it.idx() != bf_final(1))
                            bf_final(2) = v_it.idx();
;                   }
                    bco_final(2) = 0;
                }

                bco_final = bco_final / bco_final.sum();
            }

            // feature edge?
            if (efeature_ && efeature_[e])
            {
                std::cout << "insert vertex feature vertex" << std::endl;
                auto h = mesh.insert_vertex(e, epoint[e]);
                auto v = mesh.to_vertex(h);
                auto e0 = mesh.edge(h);
                auto e1 = mesh.edge(mesh.next_halfedge(h));
                vfeature_[v] = true;
                efeature_[e0] = true;
                efeature_[e1] = true;
                bary_face_prop[mesh.to_vertex(h)] = bf_final;
                bary_coord_prop[mesh.to_vertex(h)] = bco_final;
                assert(oldvert_prop[pmp::Vertex(bary_face_prop[mesh.to_vertex(h)](2))]);
            }  

            // normal edge
            else 
            {
                auto he = mesh.insert_vertex(e, epoint[e]);

                // set new vertex barycentric coordinate
                bary_face_prop[mesh.to_vertex(he)] = bf_final;
                bary_coord_prop[mesh.to_vertex(he)] = bco_final;
                assert(oldvert_prop[pmp::Vertex(bary_face_prop[mesh.to_vertex(he)](2))]);
            }
        }

        // split faces
        pmp::Halfedge h;
        for (auto f : mesh.faces())
        {
            h = mesh.halfedge(f);
            mesh.insert_edge(h, mesh.next_halfedge(mesh.next_halfedge(h)));
            h = mesh.next_halfedge(h);
            mesh.insert_edge(h, mesh.next_halfedge(mesh.next_halfedge(h)));
            h = mesh.next_halfedge(h);
            mesh.insert_edge(h, mesh.next_halfedge(mesh.next_halfedge(h)));
        }

        // clean-up properties
        mesh.remove_vertex_property(vpoint);
        mesh.remove_edge_property(epoint);
    }

    void ssp(int data_id, pmp::SurfaceMesh& coarse, pmp::SurfaceMesh& original, pmp::SurfaceMesh& mapped, const std::vector<neuralSubdiv::RandomDecimation::DecInfo>& infos)
    {
        pmp::SurfaceMesh save;
        mapped = pmp::SurfaceMesh(coarse);
        auto bary_face_prop = mapped.add_vertex_property<Eigen::Vector3i>("v:bary_face");
        auto bary_coord_prop = mapped.add_vertex_property<Eigen::Vector3d>("v:bary_coord");
        auto oldvert_prop = mapped.add_vertex_property<bool>("v:oldvert");
        auto moved_prop = mapped.add_vertex_property<bool>("v:moved");
        int idx = 0;

        Eigen::Array3d bc;
        Eigen::Array3i bf;

        // create barycentric coordinate (base level)
        for (auto v_it : mapped.vertices())
        {
            auto he = coarse.halfedges(v_it).begin();
            bc(0) = 1; bc(1) = 0; bc(2) = 0;
            bf(0) = v_it.idx();
            bf(1) = coarse.to_vertex(*he).idx();
            bf(2) = coarse.to_vertex(coarse.next_halfedge(*he)).idx();
            bary_face_prop[v_it] = bf;
            bary_coord_prop[v_it] = bc;
            oldvert_prop[v_it] = true;
            moved_prop[v_it] = false;
            ++idx;
        }

        pmp::VertexProperty<pmp::Point> points = mapped.vertex_property<pmp::Point>("v:point");
        int nbsubd = 2;
        for (int i = 0; i < nbsubd+1; ++i)
        {
            // move vertices to their position on the original mesh
            for (auto v_it : mapped.vertices())
            {
                if (moved_prop[v_it]) continue;
                bc = bary_coord_prop[v_it];
                bf = bary_face_prop[v_it];
                coarse_to_fine(bc, bf, infos);
                points[v_it] = bc(0) * original.position(pmp::Vertex(bf(0)));
                points[v_it] += bc(1) * original.position(pmp::Vertex(bf(1)));
                points[v_it] += bc(2) * original.position(pmp::Vertex(bf(2)));
                moved_prop[v_it] = true;
            }

            // copy and save mapped mesh
            save = pmp::SurfaceMesh(mapped);
            save.garbage_collection();
            std::string filename_mapped("subd");
            filename_mapped += std::to_string(i);
            filename_mapped += "_";
            filename_mapped += std::to_string(data_id);
            filename_mapped += ".obj";
            save.write(filename_mapped);

            if (i < nbsubd) 
                upsample_mid_point(mapped, coarse);
        }

        mapped.remove_vertex_property(bary_face_prop);
        mapped.remove_vertex_property(bary_coord_prop);
        mapped.remove_vertex_property(oldvert_prop);
        mapped.remove_vertex_property(moved_prop);
    }

    void coarse_to_fine(Eigen::Array3d& bary_coord,
                        Eigen::Array3i& bary_face,
                        const std::vector<neuralSubdiv::RandomDecimation::DecInfo>& infos)
    {
        neuralSubdiv::RandomDecimation::DecInfo info;
        int v0, v1, v2, count = 0;
        Eigen::Vector2d point_uv;
        double epsilon = 1e-5;
        
        for (int i = infos.size()-1; i >= 0; i--)
        {
            info = infos[i];
            if ((bary_face.array() == info.vi.idx()).any()) // is the vertex in the onering of remaining vertex of the collapse?
            {
                assert((bary_coord.array() >= 0).all());
                assert((bary_coord.array() <= 1).all());
                assert(bary_coord.sum() - 1.0 <= epsilon);
                
                v0 = v1 = v2 = -1;
                for (int j = 0; j < info.V_map_after.rows(); ++j)
                {
                    if (info.V_map_after(j) == bary_face[0]) v0 = j;
                    if (info.V_map_after(j) == bary_face[1]) v1 = j;
                    if (info.V_map_after(j) == bary_face[2]) v2 = j;
                    // TODO: stop if all idx found ? but vmap is small so not to much gain
                }
                assert(v0 != -1 && v1 != -1 && v2 != -1);
                
                // compute point position in the UV domain
                point_uv = bary_coord(0) * info.uv_after.array().row(v0)
                         + bary_coord(1) * info.uv_after.array().row(v1)
                         + bary_coord(2) * info.uv_after.array().row(v2);

                // find its barycentric coordinates for all faces (pre collapse)
                Eigen::MatrixXd bary_uv(info.F_uv.rows(), 3);
                uv_bary_coord(point_uv, info.uv, info.F_uv, bary_uv);

                // find the valid barycentric coordinates
                int uv_face_idx = -1;
                for (int j = 0; j < bary_uv.rows(); ++j)
                {
                    if (bary_uv(j, 0) >= -epsilon && bary_uv(j, 0) <= 1.0 + epsilon
                     && bary_uv(j, 1) >= -epsilon && bary_uv(j, 1) <= 1.0 + epsilon
                     && bary_uv(j, 2) >= -epsilon && bary_uv(j, 2) <= 1.0 + epsilon)
                    {
                        uv_face_idx = j;
                        break;
                    }
                }
                if (uv_face_idx == -1)
                {
                    std::cout << "ok here's the result (uv): " << neuralSubdiv::check_triangle_quality_uv(info.uv, info.F_uv) << std::endl;
                    assert(uv_face_idx != -1);
                }
                assert(uv_face_idx != -1);

                // find corresponding vertex idx in V
                bary_coord = bary_uv.row(uv_face_idx);
                for (int j = 0; j < 3; ++j)
                    bary_face(j) = info.V_map(info.F_uv(uv_face_idx, j));

                // change sign if negative
                if ((bary_coord.array() < 0).any())
                {
                    std::cout << "bary coord before sign change : " << std::endl;
                    std::cout << bary_coord << std::endl;
                    bary_coord = bary_coord.array() - bary_coord.minCoeff();
                    std::cout << "bary coord after sign change : " << std::endl;
                    std::cout << bary_coord << std::endl;
                }
                bary_coord = bary_coord / bary_coord.sum();
                ++count;
            }
        }
    }


}