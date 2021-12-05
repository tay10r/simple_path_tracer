cmake_minimum_required(VERSION 3.8)

add_library(rt "rt.cpp" "rt.h")

target_link_libraries(rt PUBLIC ${rt_glm_target})

target_include_directories(rt PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
