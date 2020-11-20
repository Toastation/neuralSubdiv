#pragma once

#include <pmp/visualization/MeshViewer.h>

using namespace pmp;

class TestViewer : public MeshViewer
{
public:
	TestViewer(const char* title, int width, int height);

protected:
	virtual void process_imgui();
};

