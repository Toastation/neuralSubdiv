#include "utils.h"
#include "RandomDecimation.h"
#include "TestViewer.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <time.h>
#include <cmath>

#include <pmp/SurfaceMesh.h>

//using namespace neuralSubdiv;

int main(int argc, char** argv) {

	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " targetMesh.obj" << std::endl;
		return 1;
	}

	std::srand(time(NULL));

	TestViewer window("Test gendata", 800, 600);

	
	pmp::SurfaceMesh mesh;
	bool res = mesh.read(argv[1]);
	std::cout << "res : " << res << std::endl;
	std::cout << "vertex count: " << mesh.n_vertices() << std::endl;

	int tv_average = 300; // target average number of vertices
	int tv_variance = 100;
	int nb_subd = 2;
	double r = std::rand() / (double)RAND_MAX;
	int nTargetVertex = tv_average + round(tv_variance * (r - 0.5));

	neuralSubdiv::normalize_unit_box(mesh);

	//neuralSubdiv::RandomDecimation rd(mesh);
	//rd.initialize();

	window.load_mesh(argv[1]);

	// random decimation

	// upsample

	// remove redundant vertices

	return window.run();
}