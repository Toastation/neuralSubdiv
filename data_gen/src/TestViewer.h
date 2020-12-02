#pragma once

#include <pmp/visualization/MeshViewer.h>

class TestViewer : public pmp::MeshViewer
{
public:
	TestViewer(const char* title, int width, int height);

protected:
	virtual void process_imgui();
};

