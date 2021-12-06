#include "example_base.h"

#include <window_blit/window_blit.hpp>

#include <GLFW/glfw3.h>

#include <cstdlib>

#ifdef _WIN32
#include <Windows.h>
#endif

namespace {

class App final : public window_blit::AppBase
{
public:
  App(GLFWwindow* window)
    : window_blit::AppBase(window)
    , m_example(Example::Create())
  {}

  void render(float* rgb, int w, int h) override
  {
    m_example->Render(rgb, glm::ivec2(0, 0), glm::ivec2(w, h), glm::ivec2(w, h));
  }

private:
  std::unique_ptr<Example> m_example;
};

} // namespace

#ifdef _WIN32
INT
WinMain(HINSTANCE, HINSTANCE, PSTR, int)
#else
int
main()
#endif
{
  return window_blit::run_glfw_window(window_blit::AppFactory<App>());
}
