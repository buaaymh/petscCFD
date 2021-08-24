link_libraries(${PETSC_LIB})
add_executable(demo_main main.cpp)
set_target_properties(demo_main PROPERTIES OUTPUT_NAME main)
