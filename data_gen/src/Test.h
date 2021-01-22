#include <pmp/SurfaceMesh.h>

namespace neuralSubdiv
{
	bool check_F_uv(pmp::SurfaceMesh& mesh, Eigen::MatrixXi& F_map, Eigen::MatrixXi& F_uv, Eigen::MatrixXi& V_map);
	bool check_F_uv_after(pmp::SurfaceMesh& mesh, pmp::Vertex vi, Eigen::MatrixXi& F_uv_after, Eigen::MatrixXi& V_map_after);
};