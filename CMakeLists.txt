cmake_minimum_required(VERSION 3.14)
set(CMAKE_C_COMPILER "/usr/local/bin/gcc-11")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")
set(CMAKE_BUILD_TYPE "Release")
project(buaaCFD)
set(CMAKE_CXX_STANDARD 17)

add_definitions(-w)

# GoogleTest related settings
option(${PROJECT_NAME}_SUBMODULE_GOOGLETEST "Add GoogleTest as a git submodule." "ON")
# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
if (${PROJECT_NAME}_SUBMODULE_GOOGLETEST)
  add_subdirectory(GoogleTest EXCLUDE_FROM_ALL)
else (${PROJECT_NAME}_SUBMODULE_GOOGLETEST)
  # Download and unpack `googletest` at configure time
  configure_file(CMakeLists.txt.in googletest/download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest/download )
  if (result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif ()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest/download )
  if (result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif ()
  # Add `googletest` directly to our build.
  # This defines the `gtest` and `gtest_main`s targets.
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest/src
                   ${CMAKE_CURRENT_BINARY_DIR}/googletest/build
                   EXCLUDE_FROM_ALL)
  # The `gtest`/`gtest_main` targets carry header search path
  # dependencies automatically when using CMake 2.8.11 or later.
  # Otherwise we have to add them here ourselves.
  if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
  endif ()
  # Now simply link against gtest or gtest_main as needed.
endif (${PROJECT_NAME}_SUBMODULE_GOOGLETEST)
# End of GoogleTest related settings

option(${PROJECT_NAME}_OPEN_MP "Build tests for this project." "OFF")
if (${PROJECT_NAME}_OPEN_MP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif ()

# add Eigen include
set(EIGEN_INC "${PROJECT_SOURCE_DIR}/eigen")
include_directories(${EIGEN_INC})

# add Petsc include and library
set(PETSC_DIR "/Users/yangminghao/code/petsc/petsc-3.15.0")
set(PETSC_ARCH "arch-darwin-c-debug")
set(PETSC_LIB "petsc" "mpich" "lapack" "blas" "parmetis" "triangle")
include_directories("${PETSC_DIR}/include"
                    "${PETSC_DIR}/${PETSC_ARCH}/include")
link_directories("${PETSC_DIR}/${PETSC_ARCH}/lib")

# Additional headers that depends on ${PROJECT_SOURCE_DIR}
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")
configure_file(
  "${PROJECT_SOURCE_DIR}/data/path.hpp.in"
  "${PROJECT_BINARY_DIR}/include/data/path.hpp"
)

option(${PROJECT_NAME}_BUILD_TESTS "Build tests for this project." "ON")
if (${PROJECT_NAME}_BUILD_TESTS)
  include(CTest)
  add_subdirectory(test)
endif ()

option(${PROJECT_NAME}_BUILD_DEMOS "Build demos for this project." "ON")
if (${PROJECT_NAME}_BUILD_DEMOS)
  add_subdirectory(demo)
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/result")
endif (${PROJECT_NAME}_BUILD_DEMOS)