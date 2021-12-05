#pragma once

#include <glm/fwd.hpp>

#include <memory>

namespace RT {

enum class TextureTarget
{
	Albedo,
	Emission
};

enum class BufferTarget
{
	/// <summary>
	/// This bind target is used when a buffer contains data for vertex attributes. 
	/// For example, if a buffer contains the position and normals of a triangle, then
	/// it is bound to this target when it is drawn.
	/// </summary>
	Attrib,
	/// <summary>
	/// This bind target is used when a buffer contains indices to vertex attributes.
	/// It is useful when the same vertex attributes are used for more than one vertex.
	/// </summary>
	Indices
};

enum class PrimitiveKind
{
	Triangles,
	Spheres
};

class Texture
{
public:
  virtual ~Texture() = default;

  virtual void Resize(int w, int h) = 0;
};

class Buffer
{
public:
  virtual ~Buffer() = default;

  virtual void Resize(int size) = 0;
};

/// <summary>
/// A geometry object is created by a draw call. Since ray tracing requires the entire scene
/// to be rendered, drawing primitives will generate a geometry object and BVH instead of rasterizing
/// the primitives to the framebuffer. Since constructing a BVH can be an expensive process, geoemtry
/// objects can be reused from frame to frame.
/// </summary>
class Geometry
{
public:
  virtual ~Geometry() = default;
};

class DebugCallback
{
public:
  virtual ~DebugCallback() = default;

  virtual void OnError(const char* msg) = 0;
};

class Context
{
public:
  static std::unique_ptr<Context> Create();

  virtual ~Context() = default;

  virtual void AddDebugCallback(DebugCallback*) = 0;

  virtual std::shared_ptr<Buffer> GenBuffer() = 0;

  virtual std::shared_ptr<Texture> GenTexture() = 0;

  virtual void Bind(TextureTarget target, std::shared_ptr<Texture> texture) = 0;

  virtual void Bind(BufferTarget target, std::shared_ptr<Buffer> buffer) = 0;

  /// <summary>
  /// Generates a geometry object and BVH from a description of primitives.
  /// </summary>
  /// <param name="primitiveKind">The kind of primitive being drawn.</param>
  /// <param name="first">The index of the first vertex in the attribute buffer.</param>
  /// <param name="count">The number of vertices in the attribute buffer to use.</param>
  /// <returns>
  /// A geometry object resulting from the draw operation.
  /// </returns>
  virtual std::shared_ptr<Geometry> Draw(PrimitiveKind primitiveKind, int first, int count) = 0;

  /// <summary>
  /// Causes the geometry object to be rendered in the next frame.
  /// </summary>
  /// <param name="geometry">The geometry to include in the next frame.</param>
  virtual void Attach(const std::shared_ptr<Geometry> &geometry) = 0;
};

} // namespace RT
