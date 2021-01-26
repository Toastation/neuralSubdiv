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
	//bool res = mesh.read(argv[1]);
	bool res = mesh.read("C:/Users/Melvin.Melvin-PC/work/neuralSubdiv/base_implementation/data_meshes/objs_original/bob.obj");
	std::cout << "Loading mesh... " << res << std::endl;
	std::cout << "Vertex count: " << mesh.n_vertices() << std::endl;

	int nb_generated_data = 5;
	int tv_average = 100; // target average number of vertices
	int tv_variance = 10; 
	int nb_subd = 2;
	int nb_target_vertex;
	double r;

	neuralSubdiv::normalize_unit_box(mesh);

	pmp::SurfaceMesh mapped;
	for (int i = 0; i < 1; ++i)
	{
		r = std::rand() / (double)RAND_MAX;
		nb_target_vertex = tv_average + round(tv_variance * (r - 0.5));
		std::cout << "Target nb of vertices: " << nb_target_vertex << std::endl;

		pmp::SurfaceMesh mesh_copy = pmp::SurfaceMesh(mesh);
		neuralSubdiv::RandomDecimation rd(mesh_copy);

		for (auto vit : mesh_copy.vertices())
			assert(mesh_copy.is_manifold(vit));
		
		rd.initialize();
		rd.simplify(400);

		/*for (int i = rd.get_dec_infos().size() - 1; i >= 0; i--)
		{
			auto info = rd.get_dec_infos()[i];
			std::cout << "vi: " << info.vi << std::endl;
		}*/
		for (auto vit : mesh_copy.vertices())
			assert(mesh_copy.is_manifold(vit));

		neuralSubdiv::ssp(i, mesh_copy, mesh, mapped, rd.get_dec_infos());
		mesh_copy.garbage_collection();
		
		std::cout << "vertices after: " << mesh_copy.n_vertices() << std::endl;
		
		//neuralSubdiv::normalize_unit_box(mapped);

		std::string filename("copy");
		filename += std::to_string(i);
		filename += ".obj";
		mesh_copy.write(filename);
		std::cout << "Data " << i << "generated as " << filename << std::endl;
	}

	neuralSubdiv::normalize_unit_box(mesh);
	std::string filename_mapped("copyoriginal");
	filename_mapped += ".obj";
	mesh.write(filename_mapped);


	return 0;
	//window.load_mesh(argv[1]);
	//return window.run();
}