add_executable(test_geometry_edge edge.cpp)
set_target_properties(test_geometry_edge PROPERTIES OUTPUT_NAME edge)
add_test(NAME TestGeometryEdge COMMAND edge)

add_executable(test_geometry_cell cell.cpp)
set_target_properties(test_geometry_cell PROPERTIES OUTPUT_NAME cell)
add_test(NAME TestGeometryCell COMMAND cell)

add_executable(test_geometry_mesh mesh.cpp)
set_target_properties(test_geometry_mesh PROPERTIES OUTPUT_NAME mesh)
add_test(NAME TestGeometryMesh COMMAND mesh)