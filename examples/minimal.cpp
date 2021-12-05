#include "example_base.h"

#include <rt.h>

namespace {

class MinimalExample final : public Example
{
public:
	MinimalExample()
	{
		auto vertexBuffer = m_context->GenBuffer();

		vertexBuffer->Resize(2);

		vertexBuffer->Write();

		m_context->Bind(RT::BufferTarget::Attrib, vertexBuffer);

		auto sphereGeom = m_context->Draw(RT::PrimitiveKind::Spheres, 0, 2);

		m_context->Bind(RT::BufferTarget::Attrib, nullptr);

		//m_context->Attach(sphereGeom);
	}

	void Render(float* rgb, glm::ivec2 min, glm::ivec2 max, glm::ivec2 frameSize) override {

          const glm::ivec2 area = max - min;

          const int pixelCount = area.x * area.y;

		  for (int i = 0; i < pixelCount; i++) {

            const int x = min.x + (i % area.x);
            const int y = min.y + (i / area.x);

			const float u = (x + 0.5f) / area.x;
			const float v = (y + 0.5f) / area.y;

			rgb[(i * 3) + 0] = u;
			rgb[(i * 3) + 1] = v;
			rgb[(i * 3) + 2] = 1;
		  }
	}

private:
  std::unique_ptr<RT::Context> m_context{ RT::Context::Create() };
};

}

std::unique_ptr<Example>
Example::Create()
{
  return std::unique_ptr<Example>(new MinimalExample());
}
