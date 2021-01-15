#include "RandomDecimation.h"

#include <vector>

#include <pmp/SurfaceMesh.h>

namespace neuralSubdiv {

	void upsample_mid_point(pmp::SurfaceMesh& mesh);

	void ssp(pmp::SurfaceMesh& mesh, neuralSubdiv::RandomDecimation::DecInfo dec_info);

	void coarse_to_fine(pmp::SurfaceMesh& coarse, std::vector<neuralSubdiv::RandomDecimation::DecInfo>& infos);
}