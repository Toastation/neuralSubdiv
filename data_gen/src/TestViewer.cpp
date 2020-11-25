#include "TestViewer.h"
#include "RandomDecimation.h"

#include <imgui.h>
#include <pmp/algorithms/SurfaceSimplification.h>

using namespace pmp;

TestViewer::TestViewer(const char* title, int width, int height) : MeshViewer(title, width, height) {
	set_draw_mode("Hidden line");
	crease_angle_ = 0.0;
}

void TestViewer::process_imgui()
{
    MeshViewer::process_imgui();

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
            neuralSubdiv::RandomDecimation rd(mesh_);
            rd.initialize();
            rd.simplify(mesh_.n_vertices() * 0.01 * target_percentage);
            update_mesh();
            std::cout << "nb vertices after decimation: " << mesh_.n_vertices() << std::endl;
        }
    }
}