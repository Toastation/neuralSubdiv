#include "Utils.h"

#include "igl/slice.h"
#include "igl/remove_unreferenced.h"
#include "igl/lscm.h"

#include "pmp/algorithms/SurfaceNormals.h"

#include <vector>
#include <set>

namespace neuralSubdiv {

	void normalize_unit_box(pmp::SurfaceMesh& meshIn) 
	{
		pmp::BoundingBox box = meshIn.bounds();
		pmp::Point min = box.min();
		double maxComponent = -10e9;
		for (auto vit : meshIn.vertices()) {
			meshIn.position(vit) -= min;
			for (int i = 0; i < 3; i++)
				if (meshIn.position(vit)[i] > maxComponent)
					maxComponent = meshIn.position(vit)[i];
		}
		for (auto vit : meshIn.vertices()) {
			meshIn.position(vit) /= maxComponent;
		}
	}


	void flatten_one_ring(pmp::SurfaceMesh& meshIn, pmp::Edge& e)
	{
		pmp::Vertex vi = meshIn.vertex(e, 0);
		pmp::Vertex vj = meshIn.vertex(e, 1);
		Eigen::MatrixXd V(meshIn.n_vertices(), 3);
		Eigen::MatrixXi F(meshIn.n_faces(), 3);
		Eigen::MatrixXi b;
		Eigen::MatrixXd bc;
		
		//auto tex = meshIn.vertex_property<pmp::TexCoord>("v:tex");

		neuralSubdiv::positions_to_matrix(meshIn.positions(), V); // TODO: should use Eigen::Map instead
		neuralSubdiv::faces_to_matrix(meshIn, F);

		// find one ring faces (without duplicates, easier with the pmp halfedge structure)
		std::vector<int> one_ring_faces(meshIn.valence(vi) + meshIn.valence(vj));
		auto it = one_ring_faces.begin();
		for (auto face : meshIn.faces(vi))
		{
			if (face.is_valid() && std::find(one_ring_faces.begin(), it, face.idx()) == it)
			{
				(*it) = face.idx();
				it++;
			}
		}
		for (auto face : meshIn.faces(vj))
		{
			if (face.is_valid() && std::find(one_ring_faces.begin(), it, face.idx()) == it)
			{
				(*it) = face.idx();
				it++;
			}
		}
		one_ring_faces.resize(std::distance(one_ring_faces.begin(), it));
		
		// get faces that are in the one-ring
		Eigen::MatrixXi F_onering_idx = Eigen::Map<Eigen::MatrixXi>(one_ring_faces.data(), one_ring_faces.size(), 1);
		Eigen::MatrixXi F_onering;
		igl::slice(F, F_onering_idx, 1, F_onering);

		// get vertices that are in the one-ring
		Eigen::MatrixXi F_uv, I, J;
		Eigen::MatrixXd V_uv;
		igl::remove_unreferenced(V, F_onering, V_uv, F_uv, I, J);

		// find vi and vj indices FUV
		int vi_uv = -1, vj_uv = -1;
		for (int i = 0; i < J.rows(); i++)
		{
			if (J(i) == static_cast<int>(vi.idx()))
				vi_uv = i;
			if (J(i) == static_cast<int>(vj.idx()))
				vj_uv = i;
			if (vi_uv != -1 && vj_uv != -1) break;
		}

		// boundary vertex idx in the uv reference
		Eigen::VectorXi boundary_idx(2);
		boundary_idx << vi_uv, vj_uv;

		double ij_norm = static_cast<double>(pmp::norm((meshIn.position(vi) - meshIn.position(vj))));

		// corresponding boundary constraints
		Eigen::MatrixXd boundary_constraints(2, 2);
		boundary_constraints << 0.0, 0.0,
								5.0*ij_norm, 0.0;

		// minimize conformal energy
		Eigen::MatrixXd uv;
		igl::lscm(V_uv, F_uv, boundary_idx, boundary_constraints, uv);
	}
}

