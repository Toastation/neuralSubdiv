#pragma once

#include <pmp/SurfaceMesh.h>

using namespace pmp;


namespace neuralSubdiv {
	/**
	 * Normalize meshIn vertex coords inside the unit box
	 *
	 * TODO: pass in a meshOut instead of modifying the meshIn?
	 */
	void normalizeUnitBox(SurfaceMesh& meshIn) {
		BoundingBox box = meshIn.bounds();
		Point min = box.min();
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
		//box = meshIn.bounds();
		//std::cout << "min " << box.min() << "  |  max " << box.max() << std::endl;
	}
} // namespace neuralSubdiv