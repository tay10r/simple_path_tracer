cmake_minimum_required(VERSION 3.14.7)

include(FetchContent)
FetchContent_Declare(windowblit URL "https://github.com/tay10r/windowblit/archive/main.zip")
FetchContent_MakeAvailable(windowblit)

add_library(rt_example_base STATIC
  "example_base.h"
  "example_base.cpp")

target_link_libraries(rt_example_base PUBLIC rt window_blit)

function(add_rt_example name)

  set(target "rt_example_${name}")

  add_executable(${target} WIN32 ${ARGN})

  target_link_libraries(${target} PRIVATE rt_example_base)

  set_target_properties(${target}
    PROPERTIES
      OUTPUT_NAME ${name}
      RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

endfunction(add_rt_example name)

add_rt_example(minimal minimal.cpp)
