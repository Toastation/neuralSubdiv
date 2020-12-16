#include "TestViewer.h"
#include "RandomDecimation.h"
#include "Utils.h"

#include <imgui.h>
#include <pmp/algorithms/SurfaceSimplification.h>
#include <pmp/algorithms/SurfaceNormals.h>

TestViewer::TestViewer(const char* title, int width, int height) : pmp::MeshViewer(title, width, height) {
	set_draw_mode("Hidden line");
	crease_angle_ = 0.0;
}

void TestViewer::process_imgui()
{
    pmp::MeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Decimation", ImGuiTreeNodeFlags_DefaultOpen))
    {
        static int target_percentage = 50;
        ImGui::PushItemWidth(100);
        ImGui::SliderInt("Percentage", &target_percentage, 1, 99);
        ImGui::PopItemWidth();

        static int normal_deviation = 180;
        ImGui::PushItemWidth(100);
        ImGui::SliderInt("Normal Deviation", &normal_deviation, 1, 180);
        ImGui::PopItemWidth();

        static int aspect_ratio = 10;
        ImGui::PushItemWidth(100);
        ImGui::SliderInt("Aspect Ratio", &aspect_ratio, 1, 10);
        ImGui::PopItemWidth();

        if (ImGui::Button("Decimate it!"))
        {
            /*neuralSubdiv::RandomDecimation rd(mesh_);
            rd.initialize();
            rd.simplify(mesh_.n_vertices() * 0.01 * target_percentage);
            update_mesh();
            std::cout << "nb vertices after decimation: " << mesh_.n_vertices() << std::endl;*/
            Eigen::MatrixXi F(mesh_.n_faces(), 3);
            neuralSubdiv::faces_to_matrix(mesh_, F);
            
            pmp::Edge edge(403);
            pmp::Vertex vi = mesh_.vertex(edge, 0); pmp::Vertex vj = mesh_.vertex(edge, 1);
            Eigen::Vector2i boundary_idx;
            Eigen::MatrixXi F_uv, F_onering, V_map;
            Eigen::MatrixXd uv, boundary_constraints;
            Eigen::ArrayXi F_map;
            neuralSubdiv::flatten_one_ring(mesh_, edge, F, 
                                           uv, F_uv, F_onering, 
                                           boundary_idx, boundary_constraints, 
                                           V_map, F_map);

            std::cout << "flatten onering ok" << std::endl;

            neuralSubdiv::check_lscm_self_folding(uv, F_uv, boundary_idx);

            std::cout << "check lscm ok" << std::endl;

            neuralSubdiv::check_link_condition(mesh_, pmp::Edge(30));
            
            std::cout << "check linked condition ok" << std::endl;

            std::vector<pmp::Normal> face_normals(F_onering.rows());
            // compute normals for face onering (normalized)
            for (int i = 0; i < F_onering.rows(); ++i)
                face_normals[i] = pmp::SurfaceNormals::compute_face_normal(mesh_, pmp::Face(F_onering(i)));
            
            // replace occurences of vj with vi
            F = (F.array() == vj.idx()).select(vi.idx(), F);
            // detect and check faces that need to be deleted
            pmp::Face del_face[2];
            int del_face_count = 0;
            for (int idx : F_map)
            {
                if (F(idx, 0) == F(idx, 1) || F(idx, 0) == F(idx, 2) || F(idx, 1) == F(idx, 2))
                {
                    if (del_face_count > 1)
                    {
                        std::cerr << "Error, #faces to delete should be 2" << std::endl;
                        break; // return 
                    }
                    del_face[del_face_count] = pmp::Face(idx);
                    del_face_count++;
                }
            }

            update_mesh();
        }
    }
}