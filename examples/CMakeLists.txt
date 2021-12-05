cmake_minimum_required(VERSION 3.14.7)

include(FetchContent)
FetchContent_Declare("glfw" URL "https://github.com/glfw/glfw/archive/master.zip")
FetchContent_MakeAvailable("glfw")

add_library(rt_example_base STATIC
	"example_base.h"
	"example_base.cpp")

target_link_libraries(rt_example_base PUBLIC rt glfw)

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
