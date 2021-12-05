#include "rt.h"

#include <glm/glm.hpp>

#include <vector>

namespace RT {

namespace {

class BufferImpl final : public Buffer
{
public:
  void Resize(int size) override { m_data.reset(new Byte[size]); }

private:
  using Byte = unsigned char;

  std::unique_ptr<Byte[]> m_data;
};

class TextureImpl final : public Texture
{
public:
  void Resize(int w, int h) override
  {
    m_width = w;
    m_height = h;
    m_color.reset(new glm::vec3[w * h]);
  }

private:
  int m_width = 0;
  int m_height = 0;
  std::unique_ptr<glm::vec3[]> m_color;
};

class TriangleGeometry;
class SphereGeometry;

class GeometryVisitor
{
public:
  virtual ~GeometryVisitor() = default;
  virtual void Visit(const TriangleGeometry&) = 0;
  virtual void Visit(const SphereGeometry&) = 0;
};

class GeometryImpl : public Geometry
{
public:
  virtual ~GeometryImpl() = default;

  virtual void Accept(GeometryVisitor& visitor) const = 0;
};

class SphereGeometry final : public GeometryImpl
{
public:
  void Accept(GeometryVisitor& visitor) const override { visitor.Visit(*this);  }
};

class TriangleGeometry final : public GeometryImpl
{
public:
  void Accept(GeometryVisitor& visitor) const override { visitor.Visit(*this);  }
};

class Scene;

class GeometryResolver final : public GeometryVisitor
{
public:
  GeometryResolver(Scene* scene)
    : m_scene(scene)
  {}

  void Visit(const TriangleGeometry&) override;

  void Visit(const SphereGeometry&) override;

private:
  Scene* m_scene;
};

class Scene final
{
public:
  void Attach(const std::shared_ptr<Geometry>& geometry)
  {
    m_attachedGeometries.emplace_back(geometry);

    ResolveDerived(static_cast<const GeometryImpl*>(geometry.get()));
  }

private:
  void ResolveDerived(const GeometryImpl* geometryImpl)
  {
    GeometryResolver resolver(this);

    geometryImpl->Accept(resolver);
  }

  friend GeometryResolver;

  void Resolve(const SphereGeometry& geom)
  {
    m_sphereGeometries.emplace_back(&geom);
  }

  void Resolve(const TriangleGeometry& geom)
  {
    m_triangleGeometries.emplace_back(&geom);
  }

private:
  /// <summary>
  ///  This vector is used just so that the lifetime of each geometry object is at least
  ///  as long as the scene. Each one of these base pointers also have a derived pointer in the
  ///  scene object.
  /// </summary>
  std::vector<std::shared_ptr<Geometry>> m_attachedGeometries;

  std::vector<const SphereGeometry*> m_sphereGeometries;

  std::vector<const TriangleGeometry*> m_triangleGeometries;
};

void
GeometryResolver::Visit(const TriangleGeometry& geom)
{
  m_scene->Resolve(geom);
}

void
GeometryResolver::Visit(const SphereGeometry& geom)
{
  m_scene->Resolve(geom);
}

struct ContextState final
{
  std::shared_ptr<Texture> albedoTexture;

  std::shared_ptr<Texture> emissionTexture;

  std::shared_ptr<Buffer> attribBuffer;

  std::shared_ptr<Buffer> indexBuffer;
};

class ContextImpl final : public Context
{
public:
  void AddDebugCallback(DebugCallback* debugCallback) override
  {
    m_debugCallbacks.emplace_back(debugCallback);
  }

  std::shared_ptr<Texture> GenTexture() override
  {
    return std::shared_ptr<Texture>(new TextureImpl);
  }

  std::shared_ptr<Buffer> GenBuffer() override
  {
    return std::shared_ptr<Buffer>(new BufferImpl);
  }

  void Bind(TextureTarget target, std::shared_ptr<Texture> texture) override
  {
    switch (target) {
      case TextureTarget::Albedo:
        GetState()->albedoTexture = texture;
        break;
      case TextureTarget::Emission:
        GetState()->emissionTexture = texture;
        break;
    }
  }

  void Bind(BufferTarget target, std::shared_ptr<Buffer> buffer) override
  {
    switch (target) {
      case BufferTarget::Attrib:
        GetState()->attribBuffer = buffer;
        break;
      case BufferTarget::Indices:
        GetState()->indexBuffer = buffer;
        break;
    }
  }

  std::shared_ptr<Geometry> Draw(PrimitiveKind kind, int first, int count) override
  {
    if (!GetState()->attribBuffer) {
      Error("No attrib buffer is bound.");
      return nullptr;    
    }

    switch (kind) {
      case PrimitiveKind::Triangles:
        break;
      case PrimitiveKind::Spheres:
        break;
    }

    return nullptr;
  }

  void Attach(const std::shared_ptr<Geometry>& geometry) override
  {
    GetScene()->Attach(geometry);
  }

private:
  void Error(const char* msg) {
    for (DebugCallback* debugCallback : m_debugCallbacks)
      debugCallback->OnError(msg);
  }

  ContextState* GetState() { return &m_state; }

  Scene* GetScene()
  {
    if (!m_scene)
      m_scene.reset(new Scene);

    return m_scene.get();
  }

private:
  std::vector<DebugCallback*> m_debugCallbacks;

  std::unique_ptr<Scene> m_scene;

  ContextState m_state;
};

} // namespace

std::unique_ptr<Context>
Context::Create()
{
  return std::unique_ptr<Context>(new ContextImpl);
}

} // namespace RT
