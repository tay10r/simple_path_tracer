#include "example_base.h"

#include <rt.hpp>

namespace {

class MinimalExample final : public Example
{
public:
  MinimalExample()
  {
    auto vertexBuffer = m_context->GenBuffer();

    // clang-format off
    const float vertices[]{
      // sphere 0 (center+radius)
      -1.0f,  0.0f,  0.0f, 1.0f,
      // sphere 1
       1.0f,  0.0f, -1.0f, 1.0f,
      // sphere 2
       0.0f,  -101.0f, 0.0f, 100.0f
    };
    // clang-format on

    vertexBuffer->Resize(sizeof(vertices));

    vertexBuffer->Write(0, vertices, sizeof(vertices));

    m_context->Bind(RT::BufferTarget::Attrib, vertexBuffer);

    m_sphereGeometry = m_context->Draw(RT::PrimitiveKind::Spheres, 0, 3);

    m_context->Bind(RT::BufferTarget::Attrib, nullptr);
  }

  void Render(float* rgb, glm::ivec2 min, glm::ivec2 max, glm::ivec2 frameSize) override
  {
    m_context->Attach(m_sphereGeometry);

    m_context->SetViewport(0, 0, frameSize.x, frameSize.y);

    m_context->SetCameraPosition(0, 0, 5);

    m_context->Render(min.x, min.y, max.x, max.y, rgb);

    m_context->ClearScene();
  }

private:
  std::unique_ptr<RT::Context> m_context{ RT::Context::Create() };

  std::shared_ptr<RT::Geometry> m_sphereGeometry;
};

} // namespace

std::unique_ptr<Example>
Example::Create()
{
  return std::unique_ptr<Example>(new MinimalExample());
}
