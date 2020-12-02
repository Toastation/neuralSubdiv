#include "Utils.h"

#include "igl/lscm.h"

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

		// find one ring faces (without duplicates)
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
		float e_norm = pmp::norm((meshIn.position(vi) - meshIn.position(vj)));
	}
}

