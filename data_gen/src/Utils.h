#pragma once

//#define DEBUG_PRINT
//#define DEBUG_PRINT_MORE

#include <pmp/SurfaceMesh.h>

#include <map>

namespace neuralSubdiv {

	// Normalize meshIn vertex coords inside the unit box
	// TODO: pass in a meshOut instead of modifying the meshIn?
	void normalize_unit_box(pmp::SurfaceMesh& meshIn);

	// Flatten the one ring of the edge V[edge[0]],V[edge[1]]
	void flatten_one_ring(pmp::SurfaceMesh& meshIn, pmp::Vertex& vi, pmp::Vertex& vj,
						Eigen::MatrixXd& uv, Eigen::MatrixXi& F_uv, Eigen::MatrixXi& F_onering,
						Eigen::Vector2i& boundary_idx, Eigen::MatrixXd& boundary_constraints,
						Eigen::MatrixXi& V_map, Eigen::MatrixXi& F_map);

	void flatten_one_ring_after(pmp::SurfaceMesh& meshIn, pmp::Vertex& vi, Eigen::MatrixXd& uv, Eigen::MatrixXi& V_map,
		Eigen::MatrixXd& uv_after, Eigen::MatrixXi& f_uv, Eigen::MatrixXi& v_idx);

	// appendix C. paragraph 3
	// check that the total angle around vi and vj is 2*pi 
	bool check_lscm_self_folding(Eigen::MatrixXd& uv, Eigen::MatrixXi& f_uv, Eigen::Vector2i& boundary_idx);

	bool check_triangle_quality_uv(Eigen::MatrixXd& uv, Eigen::MatrixXi& F_uv);

	bool check_triangle_quality(pmp::SurfaceMesh& mesh, pmp::Vertex v);

	bool reconnect_faces(Eigen::MatrixXi& F, Eigen::ArrayXi& F_map, int i, int j, Eigen::Vector2i& del_faces_idx);

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

	static inline void get_onering_faces(pmp::SurfaceMesh& mesh, pmp::Vertex& i, pmp::Vertex& j, Eigen::MatrixXi& onering_vertices_idx, Eigen::MatrixXi& onering_faces_map, Eigen::MatrixXi& onering_faces_uv) {
		std::vector<int> one_ring_faces(mesh.valence(i) + mesh.valence(j));
		auto it = one_ring_faces.begin();
		for (auto face : mesh.faces(i))
		{
			if (face.is_valid() && std::find(one_ring_faces.begin(), it, face.idx()) == it)
			{
				(*it) = face.idx();
				it++;
			}
		}
		for (auto face : mesh.faces(j))
		{
			if (face.is_valid() && std::find(one_ring_faces.begin(), it, face.idx()) == it)
			{
				(*it) = face.idx();
				it++;
			}
		}
		one_ring_faces.resize(std::distance(one_ring_faces.begin(), it));
		onering_faces_map = Eigen::Map<Eigen::ArrayXi>(one_ring_faces.data(), one_ring_faces.size()); // will the std container data be removed at the end of this func?

		// TODO: optimize? (map N numbers between 0 and |V|-1 to 0 and N-1)
		std::map<int, int> v_map;
		for (int ii = 0; ii < onering_vertices_idx.rows(); ++ii)
			v_map[onering_vertices_idx(ii)] = ii;

		onering_faces_uv = Eigen::MatrixXi(one_ring_faces.size(), 3);
		int jj = 0;
		for (int ii = 0; ii < onering_faces_map.rows(); ++ii)
		{
			jj = 0;
			for (auto vertex : mesh.vertices(pmp::Face(onering_faces_map(ii))))
			{
				onering_faces_uv(ii, jj) = v_map[vertex.idx()];
				++jj;
			}
		}
	}

	static inline void get_onering_vertices(pmp::SurfaceMesh& mesh, pmp::Vertex& i, pmp::Vertex& j, Eigen::MatrixXi& onering_vertices_idx, Eigen::MatrixXd& onering_vertices, Eigen::Vector2i& boundary_idx) {
		std::vector<int> or_idx(mesh.valence(i) + mesh.valence(j) - 2);
		auto it = or_idx.begin(); int idx = 0;
		for (auto vertex : mesh.vertices(i))
		{
			if (vertex.is_valid() && std::find(or_idx.begin(), it, vertex.idx()) == it)
			{
				if (vertex != j)
				{
					(*it) = vertex.idx();
					onering_vertices.row(idx) = static_cast<Eigen::Vector3d>(mesh.position(vertex));
					++idx;
					++it;
				}
			}
		}
		for (auto vertex : mesh.vertices(j))
		{
			if (vertex.is_valid() && std::find(or_idx.begin(), it, vertex.idx()) == it)
			{
				if (vertex != i) 
				{
					(*it) = vertex.idx();
					onering_vertices.row(idx) = static_cast<Eigen::Vector3d>(mesh.position(vertex));
					++idx;
					++it;
				}
			}
		}
		// place i,j last so we can easily ignore them later
		(*it) = i.idx();
		boundary_idx(0) = idx;
		onering_vertices.row(idx) = static_cast<Eigen::Vector3d>(mesh.position(i)); ++idx; ++it;
		(*it) = j.idx();
		boundary_idx(1) = idx; 
		onering_vertices.row(idx) = static_cast<Eigen::Vector3d>(mesh.position(j)); ++idx; ++it;
		onering_vertices_idx = Eigen::Map<Eigen::ArrayXi>(or_idx.data(), or_idx.size());
	}

} // namespace neuralSubdiv