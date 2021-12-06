cmake_minimum_required(VERSION 3.14.7)

include(FetchContent)
FetchContent_Declare(bvh URL "https://github.com/madmann91/bvh/archive/master.zip")
FetchContent_MakeAvailable(bvh)

find_package(OpenMP REQUIRED)

add_library(rt "rt.cpp" "rt.hpp")

target_link_libraries(rt PRIVATE ${rt_glm_target} bvh OpenMP::OpenMP_CXX)

target_include_directories(rt PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_features(rt PRIVATE cxx_std_17)

if(NOT MSVC)
  target_compile_options(rt PRIVATE -Wall -Wextra -Werror -Wfatal-errors -Wno-unknown-pragmas)
endif(NOT MSVC)
