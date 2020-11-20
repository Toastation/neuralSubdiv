#include "utils.h"
#include "TestViewer.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <time.h>
#include <cmath>

#include <pmp/SurfaceMesh.h>

using namespace pmp;
//using namespace neuralSubdiv;

int main(int argc, char** argv) {

	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " targetMesh.obj" << std::endl;
		return 1;
	}

	std::srand(time(NULL));

	TestViewer window("Test gendata", 800, 600);

	
	SurfaceMesh mesh;
	bool res = mesh.read(argv[1]);
	std::cout << "res : " << res << std::endl;
	std::cout << "vertex count: " << mesh.n_vertices() << std::endl;

	int tvAverage = 300; // target average number of vertices
	int tvVariance = 100;
	int nbSubd = 2;
	double r = std::rand() / (double)RAND_MAX;
	int nTargetVertex = tvAverage + round(tvVariance * (r - 0.5));

	neuralSubdiv::normalizeUnitBox(mesh);

	window.load_mesh(argv[1]);

	// random decimation

	// upsample

	// remove redundant vertices

	return window.run();
}