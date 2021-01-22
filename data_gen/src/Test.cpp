#include "Test.h"
namespace neuralSubdiv
{
    bool check_F_uv(pmp::SurfaceMesh& mesh, Eigen::MatrixXi& F_map, Eigen::MatrixXi& F_uv, Eigen::MatrixXi& V_map)
    {
        if (F_uv.rows() != F_map.rows()) return false;
        int i = 0;
        for (int j = 0; j < F_uv.rows(); ++j)
        {
            i = 0;
            for (auto vv : mesh.vertices(pmp::Face(F_map(j))))
            {
                if (V_map(F_uv(j, i)) != vv.idx()) return false;
                i++;
            }
        }
        return true;
    }

    bool check_F_uv_after(pmp::SurfaceMesh& mesh, pmp::Vertex vi, Eigen::MatrixXi& F_uv_after, Eigen::MatrixXi& V_map_after)
    {
        std::vector<int> f_onering;
        for (auto f : mesh.faces(vi))
            f_onering.push_back(f.idx());
        if (f_onering.size() != F_uv_after.rows())
        {
            std::cerr << "F_onering and F_uv_after size mismatch" << std::endl;
            return false;
        }
        int i = 0;
        for (int j = 0; j < F_uv_after.rows(); ++j)
        {
            i = 0;
            for (auto vv : mesh.vertices(pmp::Face(f_onering[j])))
            {
                if (V_map_after(F_uv_after(j, i)) != vv.idx()) return false;
                i++;
            }
        }
        return true;
    }
};