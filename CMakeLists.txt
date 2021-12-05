cmake_minimum_required(VERSION 3.14.7)

project("rt")

option(RT_FETCH_GLM "Whether or not to automatically download GLM." ON)

if(RT_FETCH_GLM)
	include(FetchContent)
	FetchContent_Declare(glm URL "https://github.com/g-truc/glm/archive/refs/tags/0.9.9.8.zip")
	FetchContent_MakeAvailable(glm)
endif(RT_FETCH_GLM)

if(NOT TARGET glm AND NOT TARGET glm::glm)
	find_package(glm REQUIRED)
endif(NOT TARGET glm AND NOT TARGET glm::glm)

if(TARGET glm::glm)
	set(rt_glm_target glm::glm)
endif(TARGET glm::glm)

if(TARGET glm)
	set(rt_glm_target glm)
endif(TARGET glm)

if(NOT DEFINED rt_glm_target)
	message(FATAL_ERROR "Failed to setup GLM")
endif(NOT DEFINED rt_glm_target)

# Include sub-projects.
add_subdirectory("rt")
add_subdirectory("examples")
