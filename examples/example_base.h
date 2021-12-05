#pragma once

#include <glm/glm.hpp>

#include <memory>

class Example
{
public:
	static std::unique_ptr<Example> Create();

	virtual ~Example() = default;

	virtual void Render(float* rgb, glm::ivec2 min, glm::ivec2 max, glm::ivec2 frameSize) = 0;
};
