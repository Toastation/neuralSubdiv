#include "Upsample.h"
#include "Utils.h"

#include <igl/find.h>

namespace neuralSubdiv {
    // based on pmp::SurfaceSubdivision::loop
    void upsample_mid_point(pmp::SurfaceMesh& mesh)
    {
        pmp::VertexProperty<pmp::Point>points_ = mesh.vertex_property<pmp::Point>("v:point");
        pmp::VertexProperty<bool> vfeature_ = mesh.get_vertex_property<bool>("v:feature");
        pmp::EdgeProperty<bool> efeature_ = mesh.get_edge_property<bool>("e:feature");

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

        for (auto e : mesh.edges())
            epoint[e] = (points_[mesh.vertex(e, 0)] + points_[mesh.vertex(e, 1)]) * pmp::Scalar(0.5);

        // insert new vertices on edges
        for (auto e : mesh.edges())
        {
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
            }

            // normal edge
            else
                mesh.insert_vertex(e, epoint[e]);
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

    void ssp(pmp::SurfaceMesh& coarse, pmp::SurfaceMesh& original, pmp::SurfaceMesh& mapped, const std::vector<neuralSubdiv::RandomDecimation::DecInfo>& infos)
    {
        Eigen::SparseMatrix<double> S(coarse.n_vertices(), coarse.n_vertices());
        S.setIdentity();

        Eigen::MatrixXd bary_coord(coarse.n_vertices(), 3);   // |V|x3 
        Eigen::MatrixXi bary_face(coarse.n_vertices(), 3);    // |V|x3 

        // NEED REWORK FOR UPSAMPLING LEVELS
        int idx = 0, n_idx = 0;
        pmp::Face face;
        for (auto v_it : coarse.vertices())
        {
            face = *(coarse.faces(v_it).begin());
            n_idx = 0;
            for (auto vn_it : coarse.vertices(face))
            {
                assert(!coarse.is_deleted(vn_it));
                assert(coarse.is_valid(vn_it));
                bary_face(idx, n_idx) = vn_it.idx();
                ++n_idx;
            }
            for (int i = 0; i < 3; ++i)
            {
                if (bary_face(idx, i) == v_it.idx())
                    bary_coord(idx, i) = 1;
                else
                    bary_coord(idx, i) = 0;
            }
            ++idx;
        }
        mapped = pmp::SurfaceMesh(coarse);
        pmp::VertexProperty<pmp::Point> points = mapped.vertex_property<pmp::Point>("v:point");
        Eigen::Array3d bc/*(bary_coord(0, 0), bary_coord(0, 1), bary_coord(0, 2))*/;
        Eigen::Array3i bf/*(bary_face(0, 0), bary_face(0, 1), bary_face(0, 2))*/;
        
        idx = 0;
        for (auto v_it : coarse.vertices())
        {
            bc = bary_coord.row(idx);
            bf = bary_face.row(idx);
            coarse_to_fine(bc, bf, infos);
            std::cout << "MOVING" << std::endl;
            points[v_it] += bc(0) * original.position(pmp::Vertex(bf(0)));
            points[v_it] += bc(1) * original.position(pmp::Vertex(bf(1)));
            points[v_it] += bc(2) * original.position(pmp::Vertex(bf(2)));
            idx++;
        }
    }

    void coarse_to_fine(Eigen::Array3d& bary_coord,
                        Eigen::Array3i& bary_face,
                        const std::vector<neuralSubdiv::RandomDecimation::DecInfo>& infos)
    {
        neuralSubdiv::RandomDecimation::DecInfo info;
        int v0, v1, v2, count = 0;
        Eigen::Vector2d point_uv;
        
        for (int i = infos.size()-1; i >= 0; i--)
        {
            info = infos[i];
            //std::cout << "vi=" << info.vi << std::endl;
            if ((bary_face.array() == info.vi.idx()).any()) // is the vertex in the onering of remaining vertex of the collapse?
            {
                //std::cout << "bary_face: " << std::endl;
                //std::cout << bary_face << std::endl;
                //std::cout << "vmap after" << std::endl;
                //std::cout << info.V_map_after << std::endl;
                assert((bary_coord.array() >= 0).all());
                assert((bary_coord.array() <= 1).all());
                assert(bary_coord.sum() - 1.0 <= std::numeric_limits<double>::epsilon());
                
                v0 = v1 = v2 = -1;
                for (int j = 0; j < info.V_map_after.rows(); ++j)
                {
                    if (info.V_map_after(j) == bary_face[0]) v0 = j;
                    if (info.V_map_after(j) == bary_face[1]) v1 = j;
                    if (info.V_map_after(j) == bary_face[2]) v2 = j;
                    // TODO: stop if all idx found ? but vmap is small so not to much gain
                }
                assert(v0 != -1 && v1 != -1 && v2 != -1);
                //if (v0 == -1 || v1 == -1 || v2 == -1) {
                //    continue;
                //}
                
                // compute point position in the UV domain
                point_uv = bary_coord(0) * info.uv_after.array().row(v0)
                         + bary_coord(1) * info.uv_after.array().row(v1)
                         + bary_coord(2) * info.uv_after.array().row(v2);

                // find its barycentric coordinates for all faces (pre collapse)
                Eigen::MatrixXd bary_uv(info.F_uv.rows(), 3);
                uv_bary_coord(point_uv, info.uv, info.F_uv, bary_uv);

                // find the valid barycentric coordinates
                int uv_face_idx = -1;
                //std::cout << "point uv" << std::endl;
                //std::cout << point_uv << std::endl;

                //std::cout << "uv" << std::endl;
                //std::cout << info.uv << std::endl;

                //std::cout << "fuv" << std::endl;
                //std::cout << info.F_uv << std::endl;

                //std::cout << "bary uv" << std::endl;
                //std::cout << bary_uv << std::endl;
                for (int j = 0; j < bary_uv.rows(); ++j)
                {
                    if (bary_uv(j, 0) >= -std::numeric_limits<double>::epsilon() && bary_uv(j, 0) <= 1.0 + std::numeric_limits<double>::epsilon()
                     && bary_uv(j, 1) >= -std::numeric_limits<double>::epsilon() && bary_uv(j, 1) <= 1.0 + std::numeric_limits<double>::epsilon()
                     && bary_uv(j, 2) >= -std::numeric_limits<double>::epsilon() && bary_uv(j, 2) <= 1.0 + std::numeric_limits<double>::epsilon())
                    {
                        uv_face_idx = j;
                        //std::cout << bary_uv.row(j) << std::endl;
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

                // TODO: divide bary_coord by sum?
                ++count;
            }
        }

        std::cout << "count" << std::endl;
        std::cout << count << std::endl;
    }


}