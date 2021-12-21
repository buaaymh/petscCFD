add_executable(test_flux_riemann riemann.cpp)
set_target_properties(test_flux_riemann PROPERTIES OUTPUT_NAME riemann)
add_test(NAME TestFluxRiemann COMMAND riemann)
