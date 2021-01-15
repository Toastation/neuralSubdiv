#include "Utils.h"
#include "RandomDecimation.h"
#include "Upsample.h"
#include "TestViewer.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <string>

#include <pmp/SurfaceMesh.h>

//using namespace neuralSubdiv;

int main(int argc, char** argv) {

	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " targetMesh.obj" << std::endl;
		return 1;
	}

	std::srand(time(NULL));

	//TestViewer window("Test gendata", 800, 600);

	
	pmp::SurfaceMesh mesh;
	bool res = mesh.read(argv[1]);
	std::cout << "Loading mesh... " << res << std::endl;
	std::cout << "Vertex count: " << mesh.n_vertices() << std::endl;

	int nb_generated_data = 5;
	int tv_average = 100; // target average number of vertices
	int tv_variance = 10; 
	int nb_subd = 2;
	int nb_target_vertex;
	double r;

	neuralSubdiv::normalize_unit_box(mesh);

	for (int i = 0; i < nb_generated_data; ++i)
	{
		r = std::rand() / (double)RAND_MAX;
		nb_target_vertex = tv_average + round(tv_variance * (r - 0.5));
		std::cout << "Target nb of vertices: " << nb_target_vertex << std::endl;

		pmp::SurfaceMesh mesh_copy = pmp::SurfaceMesh(mesh);
		neuralSubdiv::RandomDecimation rd(mesh_copy);
		
		rd.initialize();
		rd.simplify(nb_target_vertex);
		
		std::string filename("copy");
		filename += std::to_string(i);
		filename += ".obj";
		mesh_copy.write(filename);
		std::cout << "Data " << i << "generated as " << filename << std::endl;
	}


	return 0;
	//window.load_mesh(argv[1]);
	//return window.run();
}