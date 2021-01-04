#include "Utils.h"

#include "igl/slice.h"
#include "igl/remove_unreferenced.h"
#include "igl/lscm.h"
#include "igl/internal_angles.h"
#include "igl/squared_edge_lengths.h"

#include "pmp/algorithms/SurfaceNormals.h"

#include <vector>
#include <map>
#include <cassert>
#include <cmath>

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


	void flatten_one_ring(pmp::SurfaceMesh& meshIn, pmp::Vertex& vi, pmp::Vertex& vj, Eigen::MatrixXi& F,
				          Eigen::MatrixXd& uv, Eigen::MatrixXi& F_uv, Eigen::MatrixXi& F_onering,
						  Eigen::Vector2i& boundary_idx, Eigen::MatrixXd& boundary_constraints,
						  Eigen::MatrixXi& V_map, Eigen::MatrixXi& F_map)
	{
		Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor > > V(meshIn.positions()[0].data(), meshIn.positions().size(), 3);
		Eigen::MatrixXd V_uv(meshIn.valence(vi) + meshIn.valence(vj) - 2, 3);
		neuralSubdiv::get_onering_vertices(meshIn, vi, vj, V_map, V_uv, boundary_idx);
		neuralSubdiv::get_onering_faces(meshIn, vi, vj, V_map, F_map, F_uv);

#ifdef DEBUG_PRINT
		std::cout << "|V| : " << V.rows() << std::endl;
		std::cout << "|F| : " << F.rows() << std::endl;
		std::cout << "n_vertices : " << meshIn.n_vertices() << std::endl;

		std::cout << "F map" << std::endl;
		std::cout << F_map << std::endl;

		std::cout << "F uv" << std::endl;
		std::cout << F_uv << std::endl;

		std::cout << "V map" << std::endl;
		std::cout << V_map << std::endl;

		std::cout << "V_uv" << std::endl;
		std::cout << V_uv << std::endl;

		std::cout << "boundary idx" << std::endl;
		std::cout << boundary_idx << std::endl;
#endif

		double ij_norm = static_cast<double>(pmp::norm((meshIn.position(vi) - meshIn.position(vj))));

		// corresponding boundary constraints
		boundary_constraints.resize(2, 2);
		boundary_constraints << 0.0, 0.0,
								5.0 * ij_norm, 0.0;

		// minimize conformal energy
		igl::lscm(V_uv, F_uv, boundary_idx, boundary_constraints, uv);
	}

	void flatten_one_ring_after(pmp::SurfaceMesh& meshIn, pmp::Vertex& vi, Eigen::MatrixXd& uv, Eigen::MatrixXi& V_map)
	{
		Eigen::MatrixXd v(meshIn.valence(vi) + 1, 3);
		Eigen::MatrixXi v_idx(meshIn.valence(vi) + 1, 1);
		Eigen::MatrixXd boundary_constraints_after(meshIn.valence(vi), 2);
		Eigen::MatrixXi boundary_constraints_idx_after(meshIn.valence(vi), 1);

		int idx = 0;
		for (auto vertex : meshIn.vertices(vi))
		{
			v.row(idx) = static_cast<Eigen::Vector3d>(meshIn.position(vertex));
			v_idx(idx) = vertex.idx();
			++idx;
		}
		v.row(idx) = static_cast<Eigen::Vector3d>(meshIn.position(vi));
		v_idx(idx) = vi.idx();

		idx = 0;
		for (int i = 0; i < boundary_constraints_idx_after.rows(); ++i)
			boundary_constraints_idx_after(i) = i;

		std::map<int, int> uv_map;
		for (int i = 0; i < V_map.rows(); ++i)
			uv_map[V_map(i)] = i;

		for (int i = 0; i < boundary_constraints_after.rows(); ++i)
			boundary_constraints_after.row(i) = uv.row(uv_map[v_idx(boundary_constraints_idx_after(i))]);

#ifdef DEBUG_PRINT
		std::cout << "boundary idx after" << std::endl;
		std::cout << boundary_constraints_idx_after << std::endl;

		std::cout << "boundary after" << std::endl;
		std::cout << boundary_constraints_after << std::endl;
#endif
	}

	bool check_lscm_self_folding(Eigen::MatrixXd& uv, Eigen::MatrixXi& f_uv, Eigen::Vector2i& boundary_idx)
	{
		Eigen::MatrixXd sl;
		Eigen::MatrixXd angles;
		igl::squared_edge_lengths(uv, f_uv, sl);
		igl::internal_angles_using_squared_edge_lengths(sl, angles);
		double sum_angle_vi = 0.0, sum_angle_vj = 0.0;
		if (angles.rows() == f_uv.rows()) return false;
		// sum angles where vi or vj appear
		for (int i = 0; i < angles.rows(); ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				if (f_uv(i, j) == boundary_idx[0])
					sum_angle_vi += angles(i, j);
				if (f_uv(i, j) == boundary_idx[1])
					sum_angle_vj += angles(i, j);
			}
		}
		return (sum_angle_vi < 2.0 * M_PI + 1e-7)
			&& (sum_angle_vj < 2.0 * M_PI + 1e-7);
	}

	bool check_link_condition(pmp::SurfaceMesh& meshIn, pmp::Edge& e)
	{
		pmp::Vertex vi = meshIn.vertex(e, 0);
		pmp::Vertex vj = meshIn.vertex(e, 1);
		std::vector<int> vi_neighbordhood;
		int vertices[2];
		int count = 0;
		// vi neighborhood vertices
		for (auto v_it : meshIn.vertices(vi))
			vi_neighbordhood.push_back(v_it.idx());
		// check that vj's has at most 2 in common with vi's
		for (auto v_it : meshIn.vertices(vj))
		{
			if (std::find(vi_neighbordhood.begin(), vi_neighbordhood.end(), v_it.idx()) != vi_neighbordhood.end())
			{
				return false; 
				vertices[count] = v_it.idx();
				++count;
			}
		}
		// check if the two common vertices are connected
		for (auto v_it : meshIn.vertices(pmp::Vertex(vertices[0])))
			if (v_it.idx() == vertices[1])
				return false;
		return true;
	}

	bool reconnect_faces(Eigen::MatrixXi& F, Eigen::ArrayXi& F_map, int i, int j, Eigen::Vector2i& del_faces_idx)
	{
		// replace occurences of vj with vi
		F = (F.array() == j).select(i, F);
		// detect and check faces that need to be deleted
		int del_face_count = 0;
		for (int idx : F_map)
		{
			if (F(idx, 0) == F(idx, 1) || F(idx, 0) == F(idx, 2) || F(idx, 1) == F(idx, 2))
			{
				if (del_face_count > 1)
				{
					std::cerr << "Error, #faces to delete should be 2" << std::endl;
					return false;
				}
				del_faces_idx[del_face_count] = idx;
				del_face_count++;
			}
		}
		return true;
	}
}

