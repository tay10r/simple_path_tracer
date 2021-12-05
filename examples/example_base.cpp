#include "example_base.h"

#include <GLFW/glfw3.h>

#include <cstdlib>

#ifdef _WIN32
#include <Windows.h>
#endif

#ifdef _WIN32
INT
WinMain(HINSTANCE, HINSTANCE, PSTR, int)
#else
int main()
#endif
{
	if (glfwInit() != GLFW_TRUE)
		return EXIT_FAILURE;

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

	GLFWwindow* window = glfwCreateWindow(640, 480, "", nullptr, nullptr);
	if (!window) {
		glfwTerminate();
		return EXIT_FAILURE;
	}

	glfwMakeContextCurrent(window);

	glfwSwapInterval(1);

	std::unique_ptr<Example> example = Example::Create();

	while (!glfwWindowShouldClose(window)) {

		glfwPollEvents();
	}

	glfwDestroyWindow(window);

	glfwTerminate();

	return EXIT_SUCCESS;
}