file(GLOB SOURCES ./*.cpp)
file(GLOB HEADERS ./*.h)

add_executable(dataGen ${SOURCES} ${HEADERS})
target_link_libraries(dataGen pmp_vis igl::core)

# if (EMSCRIPTEN)
#     set_target_properties(myviewer PROPERTIES LINK_FLAGS "--shell-file ${PROJECT_SOURCE_DIR}/external/pmp-library/src/apps/data/shell.html --preload-file ${PROJECT_SOURCE_DIR}/external/pmp-library/external/pmp-data/off/bunny.off@input.off")
# endif()