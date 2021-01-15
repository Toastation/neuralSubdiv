#include "Upsample.h"
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

    void ssp(pmp::SurfaceMesh& mesh, neuralSubdiv::RandomDecimation::DecInfo dec_info)
    {

    }

    void coarse_to_fine(pmp::SurfaceMesh& coarse, std::vector<neuralSubdiv::RandomDecimation::DecInfo>& infos)
    {

    }
}