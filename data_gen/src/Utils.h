#pragma once

#include <pmp/SurfaceMesh.h>

using namespace pmp;


namespace neuralSubdiv {

	// Normalize meshIn vertex coords inside the unit box
	// TODO: pass in a meshOut instead of modifying the meshIn?
	void normalize_unit_box(SurfaceMesh& meshIn);

	// Flatten the one ring of the edge V[edge[0]],V[edge[1]]
	void flatten_one_ring(SurfaceMesh& meshIn, Edge& edge);


} // namespace neuralSubdiv