link_libraries(${PETSC_LIB})
# add_executable(demo_periodic periodic.cpp)
# set_target_properties(demo_periodic PROPERTIES OUTPUT_NAME periodic)

# add_executable(demo_tube tube.cpp)
# set_target_properties(demo_tube PROPERTIES OUTPUT_NAME tube)

# add_executable(demo_stage stage.cpp)
# set_target_properties(demo_stage PROPERTIES OUTPUT_NAME stage)

# add_executable(demo_shock shock.cpp)
# set_target_properties(demo_shock PROPERTIES OUTPUT_NAME shock)

add_executable(demo_bubble bubble.cpp)
set_target_properties(demo_bubble PROPERTIES OUTPUT_NAME bubble)
