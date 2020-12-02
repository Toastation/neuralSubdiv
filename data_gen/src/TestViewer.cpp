#include "TestViewer.h"
#include "RandomDecimation.h"
#include "Utils.h"

#include <imgui.h>
#include <pmp/algorithms/SurfaceSimplification.h>

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
            pmp::Edge edge(403);
            Eigen::Vector2i boundary_idx;
            Eigen::MatrixXi F_uv, V_map, F_map;
            Eigen::MatrixXd uv, boundary_constraints;
            neuralSubdiv::flatten_one_ring(mesh_, edge, uv, F_uv, boundary_idx, boundary_constraints, 
                                           V_map, F_map);
            neuralSubdiv::check_lscm_self_folding(uv, F_uv, boundary_idx);
            update_mesh();
        }
    }
}