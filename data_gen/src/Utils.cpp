#include "Utils.h"

#include "igl/slice.h"
#include "igl/remove_unreferenced.h"
#include "igl/lscm.h"
#include "igl/internal_angles.h"
#include "igl/squared_edge_lengths.h"

#include "pmp/algorithms/SurfaceNormals.h"
#include "pmp/algorithms/DifferentialGeometry.h"

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


	void flatten_one_ring(pmp::SurfaceMesh& meshIn, pmp::Vertex& vi, pmp::Vertex& vj,
				          Eigen::MatrixXd& uv, Eigen::MatrixXi& F_uv, Eigen::MatrixXi& F_onering,
						  Eigen::Vector2i& boundary_idx, Eigen::MatrixXd& boundary_constraints,
						  Eigen::MatrixXi& V_map, Eigen::MatrixXi& F_map)
	{
		Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor > > V(meshIn.positions()[0].data(), meshIn.positions().size(), 3);
		Eigen::MatrixXd V_uv(meshIn.valence(vi) + meshIn.valence(vj) - 2, 3);
		neuralSubdiv::get_onering_vertices(meshIn, vi, vj, V_map, V_uv, boundary_idx);
		neuralSubdiv::get_onering_faces(meshIn, vi, vj, V_map, F_map, F_uv);

#ifdef DEBUG_PRINT_MORE
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
		assert(V_map(boundary_idx[0]) == vi.idx());
	}

	void flatten_one_ring_after(pmp::SurfaceMesh& meshIn, pmp::Vertex& vi, Eigen::MatrixXd& uv, Eigen::MatrixXi& V_map, 
		Eigen::MatrixXd& uv_after, Eigen::MatrixXi& f_uv, Eigen::MatrixXi& v_idx)
	{
		Eigen::MatrixXd v(meshIn.valence(vi) + 1, 3);
		Eigen::MatrixXd boundary_constraints_after(meshIn.valence(vi), 2);
		Eigen::MatrixXi boundary_constraints_idx_after(meshIn.valence(vi), 1);

		std::map<int, int> uv_map, v_map_rev;
		for (int i = 0; i < V_map.rows(); ++i)
			uv_map[V_map(i)] = i;

		int idx = 0;
		for (auto vertex : meshIn.vertices(vi))
		{
			v.row(idx) = static_cast<Eigen::Vector3d>(meshIn.position(vertex));
			v_idx(idx) = vertex.idx();
			v_map_rev[vertex.idx()] = idx;
			++idx;
		}
		v.row(idx) = static_cast<Eigen::Vector3d>(meshIn.position(vi));
		v_idx(idx) = vi.idx();
		v_map_rev[vi.idx()] = idx;

		int n_faces = 0;
		for (auto face : meshIn.faces(vi))
			++n_faces;

		f_uv = Eigen::MatrixXi(n_faces, 3);
		int i = 0, j = 0;
		for (auto face : meshIn.faces(vi))
		{
			j = 0;
			for (auto vertex : meshIn.vertices(face))
			{
				f_uv(i, j) = v_map_rev[vertex.idx()];
				++j;
			}
			++i;
		}

		idx = 0;
		for (int i = 0; i < boundary_constraints_idx_after.rows(); ++i)
			boundary_constraints_idx_after(i) = i;

		for (int i = 0; i < boundary_constraints_after.rows(); ++i)
			boundary_constraints_after.row(i) = uv.row(uv_map[v_idx(boundary_constraints_idx_after(i))]);

#ifdef DEBUG_PRINT_MORE
		std::cout << "v idx" << std::endl;
		std::cout << v_idx << std::endl;
		
		std::cout << "f uv" << std::endl;
		std::cout << f_uv << std::endl;

		std::cout << "boundary idx after" << std::endl;
		std::cout << boundary_constraints_idx_after << std::endl;

		std::cout << "boundary after" << std::endl;
		std::cout << boundary_constraints_after << std::endl;
#endif

		igl::lscm(v, f_uv, boundary_constraints_idx_after, boundary_constraints_after, uv_after);
		//// temp fix, WHY is it swapped ??
		for (int i = 0; i < uv_after.rows(); ++i)
			std::swap(uv_after(i, 0), uv_after(i, 1));
	}

	bool check_lscm_self_folding(Eigen::MatrixXd& uv, Eigen::MatrixXi& f_uv, Eigen::Vector2i& boundary_idx)
	{
		Eigen::MatrixXd sl;
		Eigen::MatrixXd angles;
		igl::squared_edge_lengths(uv, f_uv, sl);
		igl::internal_angles_using_squared_edge_lengths(sl, angles);
		double sum_angle_vi = 0.0, sum_angle_vj = 0.0;
		if (angles.rows() != f_uv.rows()) return false;
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

	bool check_triangle_quality_uv(Eigen::MatrixXd& uv, Eigen::MatrixXi& F_uv)
	{
		// TODO optimize and factorize
		Eigen::Vector2d p1, p2, p3;
		double el1, el2, el3, x, delta, Q;
		double cst = 4.0 * std::sqrt(3.0);
		for (int i = 0; i < F_uv.rows(); ++i)
		{
			p1 = uv.row(F_uv(i, 0));
			p2 = uv.row(F_uv(i, 1));
			p3 = uv.row(F_uv(i, 2));
			el1 = (p1 - p2).norm();
			el2 = (p2 - p3).norm();
			el3 = (p3 - p1).norm();
			x = (el1 + el2 + el3) / 2.0;
			delta = std::sqrt(x * (x - el1) * (x - el2) * (x - el3));
			Q = cst * delta / (el1 * el1 + el2 * el2 + el3 * el3);
			//std::cout << "cst:" << cst << "    denom: " << (el1 * el1 + el2 * el2 + el3 * el3) << std::endl;
			//std::cout << "Q=" << Q << "  x=" << x << "   delta=" << delta << std::endl;
			if (Q < 0.2) return false;
		}
		return true;
	}

	bool check_triangle_quality(pmp::SurfaceMesh& mesh, pmp::Vertex v)
	{
		// TODO optimize and factorize
		pmp::Point p1, p2, p3;
		double el1, el2, el3, x, delta, Q;
		double cst = 4.0 * std::sqrt(3.0);
		for (auto face : mesh.faces(v))
		{
			auto it = mesh.vertices(face);
			p1 = mesh.position(*it); ++it;
			p2 = mesh.position(*it); ++it;
			p3 = mesh.position(*it);
			el1 = dot(p1, p2);
			el2 = dot(p2, p3);
			el3 = dot(p3, p1);
			x = (el1 + el2 + el3) / 2.0;
			delta = std::sqrt(x * (x - el1) * (x - el2) * (x - el3));
			Q = cst * delta / (el1 * el1 + el2 * el2 + el3 * el3);
			//std::cout << "cst:" << cst << "    denom: " << (el1 * el1 + el2 * el2 + el3 * el3) << std::endl;
			//std::cout << "Q=" << Q << "  x=" << x << "   delta=" << delta << std::endl;
			if (Q < 0.2) return false;
		}
		return true;
	}

	/*
	l0 = sqrt(sum((UV(FUV(:, 1), :) - UV(FUV(:, 2), :)). ^ 2, 2));
	l1 = sqrt(sum((UV(FUV(:, 2), :) - UV(FUV(:, 3), :)). ^ 2, 2));
	l2 = sqrt(sum((UV(FUV(:, 3), :) - UV(FUV(:, 1), :)). ^ 2, 2));
	x = (l0 + l1 + l2) . / 2;
	delta = sqrt(x .* (x - l0) .* (x - l1) .* (x - l2));
	triQ = 4 * sqrt(3) * delta . / (l0. ^ 2 + l1. ^ 2 + l2. ^ 2);
	triQ(delFUV) = [];
	if any(triQ < triQ_threshold)
		fprintf('bad UV triangle quality\n');
	V = Vpre;
	F = Fpre;
	ECost(e) = inf;
	if stall > maxStall break; end
		continue;
	end*/

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

