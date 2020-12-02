#pragma once

#include <pmp/SurfaceMesh.h>

namespace neuralSubdiv {

	// Normalize meshIn vertex coords inside the unit box
	// TODO: pass in a meshOut instead of modifying the meshIn?
	void normalize_unit_box(pmp::SurfaceMesh& meshIn);

	// Flatten the one ring of the edge V[edge[0]],V[edge[1]]
	void flatten_one_ring(pmp::SurfaceMesh& meshIn, pmp::Edge& edge);

	// TO FIX: convert directly from meshIn.positions() to eigen matrix but seems rather tricky with Eigen::Map(?)
	// std::vector<Eigen::Matrix<float, 3, 1>>  to something like   Eigen::MatrixXf(N, 3)  
	// (it should technically also cast the scalars to double but we can take of that afterwards) 
	static inline void positions_to_matrix(std::vector<pmp::Point>& P, Eigen::MatrixXd& V) 
	{
		for (int i = 0; i < P.size(); i++)
			V.row(i) = (static_cast<Eigen::Matrix<float, 3, 1>>(P[i])).cast<double>();
	}

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