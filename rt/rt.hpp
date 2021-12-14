#pragma once

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
  /// This bind target is used when a buffer contains data for vertex
  /// attributes. For example, if a buffer contains the position and normals of
  /// a triangle, then it is bound to this target when it is drawn.
  /// </summary>
  Attrib,
  /// <summary>
  /// This bind target is used when a buffer contains indices to vertex
  /// attributes. It is useful when the same vertex attributes are used for more
  /// than one vertex.
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

  virtual void Write(int offset, const void* data, int size) = 0;
};

/// <summary>
/// A geometry object is created by a draw call. Since ray tracing requires the
/// entire scene to be rendered, drawing primitives will generate a geometry
/// object and BVH instead of rasterizing the primitives to the framebuffer.
/// Since constructing a BVH can be an expensive process, geoemtry objects can
/// be reused from frame to frame.
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
  /// <param name="first">The index of the first vertex in the attribute
  /// buffer.</param> <param name="count">The number of vertices in the
  /// attribute buffer to use.</param> <returns> A geometry object resulting
  /// from the draw operation.
  /// </returns>
  virtual std::shared_ptr<Geometry> Draw(PrimitiveKind primitiveKind, int first, int count) = 0;

  /// <summary>
  /// Causes the geometry object to be rendered in the next frame.
  /// </summary>
  /// <param name="geometry">The geometry to include in the next frame.</param>
  virtual void Attach(const std::shared_ptr<Geometry>& geometry) = 0;

  /// <summary>
  /// Renders all attached geometry.
  /// </summary>
  virtual void Render(int xMin, int yMin, int xMax, int yMax, float* rgb) = 0;

  /// <summary>
  /// Sets the transformation of screen space coordinates to NDC coordinates.
  /// </summary>
  virtual void SetViewport(int x, int y, int w, int h) = 0;

  /// <summary>
  /// Causes new random values to be drawn for each unique pixel on the next rendering operation.
  /// When multi sampling, this should be called between rendering operations on the same area.
  /// </summary>
  virtual void AdvanceRng() = 0;

  /// <summary>
  /// Constructs the top level BVH structure. Should be called once all the geometry and lighting is attached.
  /// </summary>
  virtual void CommitScene() = 0;

  /// <summary>
  /// Removes all attached geometry and lights from the scene.
  /// </summary>
  virtual void ClearScene() = 0;

  /// <summary>
  /// Sets the camera position in world space.
  /// </summary>
  virtual void SetCameraPosition(float x, float y, float z) = 0;

  /// <summary>
  /// Sets the camera rotation in radians.
  /// </summary>
  virtual void SetCameraRotation(float altitude, float azimuth) = 0;

  /// <summary>
  /// Sets the perspective transformation applies to camera rays.
  /// </summary>
  virtual void SetPerspective(float fovy, float aspect, float near, float far) = 0;
};

} // namespace RT
