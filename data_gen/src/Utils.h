#pragma once

#include <pmp/SurfaceMesh.h>

namespace neuralSubdiv {

	// Normalize meshIn vertex coords inside the unit box
	// TODO: pass in a meshOut instead of modifying the meshIn?
	void normalize_unit_box(pmp::SurfaceMesh& meshIn);

	// Flatten the one ring of the edge V[edge[0]],V[edge[1]]
	void flatten_one_ring(pmp::SurfaceMesh& meshIn, pmp::Edge& e, Eigen::MatrixXi& F,
						Eigen::MatrixXd& uv, Eigen::MatrixXi& F_uv, Eigen::MatrixXi& F_onering,
						Eigen::Vector2i& boundary_idx, Eigen::MatrixXd& boundary_constraints,
						Eigen::MatrixXi& V_map, Eigen::ArrayXi& F_map);

	// appendix C. paragraph 3
	// check that the total angle around vi and vj is 2*pi 
	bool check_lscm_self_folding(Eigen::MatrixXd& uv, Eigen::MatrixXi& f_uv, Eigen::Vector2i& boundary_idx);

	// appendix C. paragraph 4
	// check if there are non-manifold edges
	bool check_link_condition(pmp::SurfaceMesh& meshIn, pmp::Edge& e);

	// This one I don't how to use the data from pmp::SurfaceMesh 
	static inline void faces_to_matrix(pmp::SurfaceMesh& meshIn, Eigen::MatrixXi& F)
	{
		int i = 0, j = 0;
		for (pmp::Face f : meshIn.faces())
		{
			j = 0;
			for (pmp::Vertex v : meshIn.vertices(f))
			{
				F(i, j) = v.idx();
				j++;
			}
			i++;
		}
	}

} // namespace neuralSubdiv