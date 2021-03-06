#pragma once

#include <set>
#include <vector>
#include <random>

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/Heap.h"
#include "pmp/algorithms/NormalCone.h"
#include "pmp/algorithms/Quadric.h"


using namespace pmp;

namespace neuralSubdiv {



// Based on pmp::SurfaceSimplification
// Adds random decimation and succesive self-parametrization
// TODO: check if mesh is closed, check edge intersection after decimation?
class RandomDecimation {
public:
    struct DecInfo
    {
        int n_collapse;
        pmp::Vertex vi;                 // remaining vertex after the collapse
        Eigen::Vector2i boundary_idx;   // vertex idx of 1-ring
        Eigen::MatrixXd uv;             // flattened vertex position 
        Eigen::MatrixXd uv_after;       // flattened vertex position (after second collapse)
        Eigen::MatrixXi F_uv;           // faces in the uv matrix
        Eigen::MatrixXi F_uv_after;     // faces in the uv_after matrix
        Eigen::MatrixXi V_map;          // mapping from uv matrix to the mesh vertex matrix
        Eigen::MatrixXi V_map_after;
    };

    //! Construct with mesh to be simplified.
    RandomDecimation(SurfaceMesh& mesh);

    //! Destructor.
    ~RandomDecimation();

    //! Initialize with given parameters.
    void initialize(Scalar aspect_ratio = 10.0, Scalar edge_length = 0.0,
        unsigned int max_valence = 0, Scalar normal_deviation = 100.0,
        Scalar hausdorff_error = 0.0, unsigned int edge_subset_size=100, bool use_subset=true);

    //! Simplify mesh to \p n_vertices.
    void simplify(unsigned int n_vertices);

    const std::vector<neuralSubdiv::RandomDecimation::DecInfo>& get_dec_infos()
    {
        return dec_infos_;
    }

private:

    // Store data for an halfedge collapse
    struct CollapseData
    {
        CollapseData(SurfaceMesh& sm, Halfedge h);

        SurfaceMesh& mesh;

        //        vl
        //        *
        //       / \
        //      /   \
        //     / fl  \
        // v0 *------>* v1
        //     \ fr  /
        //      \   /
        //       \ /
        //        *
        //        vr
        Halfedge v0v1; // Halfedge to be collapsed
        Halfedge v1v0; // Reverse halfedge
        Vertex v0;     // Vertex to be removed
        Vertex v1;     // Remaining vertex
        Face fl;       // Left face
        Face fr;       // Right face
        Vertex vl;     // Left vertex
        Vertex vr;     // Right vertex
        Halfedge v1vl, vlv0, v0vr, vrv1;
    };

    // Heap interface
    class HeapInterface
    {
    public:
        HeapInterface(VertexProperty<float> prio, VertexProperty<int> pos)
            : prio_(prio), pos_(pos)
        {
        }

        bool less(Vertex v0, Vertex v1) { return prio_[v0] < prio_[v1]; }
        bool greater(Vertex v0, Vertex v1) { return prio_[v0] > prio_[v1]; }
        int get_heap_position(Vertex v) { return pos_[v]; }
        void set_heap_position(Vertex v, int pos) { pos_[v] = pos; }

    private:
        VertexProperty<float> prio_;
        VertexProperty<int> pos_;
    };

    typedef Heap<Vertex, HeapInterface> PriorityQueue;

    typedef std::vector<Point> Points;

    // put the vertex v in the priority queue
    void enqueue_vertex(Vertex v);

    // is collapsing the halfedge h allowed?
    bool is_collapse_legal(const CollapseData& cd);

    // what is the priority of collapsing the halfedge h
    float priority(const CollapseData& cd);

    // postprocess halfedge collapse
    void postprocess_collapse(const CollapseData& cd);

    // compute aspect ratio for face f
    Scalar aspect_ratio(Face f) const;

    // compute distance from point p to triangle f
    Scalar distance(Face f, const Point& p) const;

    // get a random subset of the mesh of size n
    void get_random_edge_subset(int n, std::vector<int>& edges_subset);

    SurfaceMesh& mesh_;

    bool initialized_;

    VertexProperty<float> vpriority_;
    VertexProperty<Halfedge> vtarget_;
    VertexProperty<int> heap_pos_;
    VertexProperty<Quadric> vquadric_;
    FaceProperty<NormalCone> normal_cone_;
    FaceProperty<Points> face_points_;

    VertexProperty<Point> vpoint_;
    FaceProperty<Point> fnormal_;
    VertexProperty<bool> vselected_;
    VertexProperty<bool> vfeature_;
    EdgeProperty<bool> efeature_;
    EdgeProperty<bool> eselected_;
    HalfedgeProperty<bool> eskip_;

    PriorityQueue* queue_;

    bool has_selection_;
    bool has_features_;
    bool use_subset_;
    Scalar normal_deviation_;
    Scalar hausdorff_error_;
    Scalar aspect_ratio_;
    Scalar edge_length_;
    unsigned int max_valence_;
    unsigned int edge_subset_size_;

    std::default_random_engine eng_;
    std::vector<neuralSubdiv::RandomDecimation::DecInfo> dec_infos_;
    std::vector<int> edges_copy;
};

} // namespace neuralSubdiv