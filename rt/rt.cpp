#include "rt.hpp"

#include <bvh/bvh.hpp>
#include <bvh/primitive_intersectors.hpp>
#include <bvh/single_ray_traverser.hpp>
#include <bvh/sphere.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/triangle.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <limits>
#include <optional>
#include <random>
#include <vector>

#include <cstring>

namespace RT {

namespace {

class BufferImpl final : public Buffer
{
public:
  const float* GetPtr() const noexcept { return reinterpret_cast<const float*>(m_data.get()); }

  int GetSize() const noexcept { return m_size; }

  void Resize(int size) override
  {
    m_data.reset(new Byte[size]);

    m_size = size;
  }

  void Write(int offset, const void* data, int size) override
  {
    if ((offset + size) > m_size)
      return;

    std::memcpy(m_data.get() + offset, data, size);
  }

private:
  using Byte = unsigned char;

  std::unique_ptr<Byte[]> m_data;

  int m_size = 0;
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
  using BvhVec3 = bvh::Vector3<float>;

  using Bvh = bvh::Bvh<float>;

  using Sphere = bvh::Sphere<float>;

  using ClosestIntersector = bvh::ClosestPrimitiveIntersector<Bvh, Sphere, false>;

  using BvhRay = bvh::Ray<float>;

  using Hit = typename ClosestIntersector::Result;

  SphereGeometry(int count)
    : m_spheres(count)
  {}

  void Accept(GeometryVisitor& visitor) const override { visitor.Visit(*this); }

  void Commit()
  {
    const auto [bboxes, centers] = bvh::compute_bounding_boxes_and_centers(m_spheres.data(), m_spheres.size());

    const auto global_bbox = bvh::compute_bounding_boxes_union(bboxes.get(), m_spheres.size());

    bvh::SweepSahBuilder<Bvh> builder(m_bvh);

    builder.build(global_bbox, bboxes.get(), centers.get(), m_spheres.size());
  }

  std::optional<Hit> FindClosestHit(const BvhRay& ray) const
  {
    ClosestIntersector primitive_intersector(m_bvh, m_spheres.data());

    bvh::SingleRayTraverser<Bvh> traverser(m_bvh);

    return traverser.traverse(ray, primitive_intersector);
  }

  void SetCenter(int index, const BvhVec3& center) { m_spheres[index].origin = center; }

  void SetRadius(int index, float radius) { m_spheres[index].radius = radius; }

  template<typename Index>
  const Sphere& operator[](Index index) const noexcept
  {
    return m_spheres[index];
  }

private:
  std::vector<Sphere> m_spheres;

  Bvh m_bvh;
};

struct Triangle final
{
  glm::vec3 p0;
  glm::vec3 e0;
  glm::vec3 e1;
  glm::vec3 n;
};

class TriangleGeometry final : public GeometryImpl
{
public:
  void Accept(GeometryVisitor& visitor) const override { visitor.Visit(*this); }
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
  using BvhVec3 = bvh::Vector<float, 3>;

  using BvhRay = bvh::Ray<float>;

  struct Hit final
  {
    PrimitiveKind kind;

    int geomIndex = -1;
    int primIndex = -1;

    float distance = std::numeric_limits<float>::infinity();

    float u = 0;
    float v = 0;

    BvhVec3 normal{ 0, 0, 0 };

    Hit() = default;

    Hit(const SphereGeometry::Hit& hit, const BvhRay& ray, const SphereGeometry& geom, int geomIndex_)
      : kind(PrimitiveKind::Spheres)
      , geomIndex(geomIndex_)
      , primIndex(hit.primitive_index)
      , distance(hit.intersection.distance())
      , u(0)
      , v(0)
      , normal(((ray.origin + (ray.direction * hit.intersection.distance())) - geom[hit.primitive_index].origin) *
               (1.0f / geom[hit.primitive_index].radius))
    {}

    bool operator>(const SphereGeometry::Hit& sphereHit) const { return distance > sphereHit.intersection.distance(); }

    constexpr operator bool() const noexcept { return distance < std::numeric_limits<float>::infinity(); }
  };

  void Attach(const std::shared_ptr<Geometry>& geometry)
  {
    m_attachedGeometries.emplace_back(geometry);

    ResolveDerived(static_cast<const GeometryImpl*>(geometry.get()));
  }

  void Commit() {}

  Hit FindClosestHit(const BvhRay& ray) const
  {
    Hit closestHit;

    for (const TriangleGeometry* geom : m_triangleGeometries) {
      (void)geom;
    }

    for (int i = 0; i < int(m_sphereGeometries.size()); i++) {

      const SphereGeometry* geom = m_sphereGeometries[i];

      const std::optional<SphereGeometry::Hit> sphereHit = geom->FindClosestHit(ray);
      if (!sphereHit)
        continue;

      if (sphereHit && (closestHit > sphereHit.value()))
        closestHit = Hit(sphereHit.value(), ray, *m_sphereGeometries[i], i);
    }

    return closestHit;
  }

private:
  void ResolveDerived(const GeometryImpl* geometryImpl)
  {
    GeometryResolver resolver(this);

    geometryImpl->Accept(resolver);
  }

  friend GeometryResolver;

  void Resolve(const SphereGeometry& geom) { m_sphereGeometries.emplace_back(&geom); }

  void Resolve(const TriangleGeometry& geom) { m_triangleGeometries.emplace_back(&geom); }

private:
  /// <summary>
  ///  This vector is used just so that the lifetime of each geometry object is
  ///  at least as long as the scene. Each one of these base pointers also have
  ///  a derived pointer in the scene object.
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

struct Viewport final
{
  int x = 0;
  int y = 0;
  int w = 1;
  int h = 1;
};

struct Perspective final
{
  float fovy = 1.0f;
  float aspect = 1.0f;
  float near = 0.0f;
  float far = std::numeric_limits<float>::infinity();
};

class ContextImpl final : public Context
{
public:
  void AddDebugCallback(DebugCallback* debugCallback) override { m_debugCallbacks.emplace_back(debugCallback); }

  std::shared_ptr<Texture> GenTexture() override { return std::shared_ptr<Texture>(new TextureImpl); }

  std::shared_ptr<Buffer> GenBuffer() override { return std::shared_ptr<Buffer>(new BufferImpl); }

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

    auto buffer = static_cast<const BufferImpl*>(GetState()->attribBuffer.get());

    auto bufferSize = buffer->GetSize();

    auto attrib = buffer->GetPtr();

    switch (kind) {
      case PrimitiveKind::Triangles:
        return DrawTriangles(first, count, attrib, bufferSize);
      case PrimitiveKind::Spheres:
        return DrawSpheres(first, count, attrib, bufferSize);
    }

    return nullptr;
  }

  void Attach(const std::shared_ptr<Geometry>& geometry) override { GetScene()->Attach(geometry); }

  void Render(int xMin, int yMin, int xMax, int yMax, float* rgb) override
  {
    const int w = xMax - xMin;
    const int h = yMax - yMin;

    const int area = w * h;

#pragma omp parallel for

    for (int i = 0; i < area; i++) {

      const int x = xMin + (i % w);
      const int y = yMin + (i / w);

      std::seed_seq seed{ x, y, m_pixelSeed };

      std::minstd_rand pixelRng(seed);

      const BvhVec3 color = RenderPixel(x, y, pixelRng);

      rgb[(i * 3) + 0] = color[0];
      rgb[(i * 3) + 1] = color[1];
      rgb[(i * 3) + 2] = color[2];
    }
  }

  void SetViewport(int x, int y, int w, int h) override { m_viewport = Viewport{ x, y, w, h }; }

  void AdvanceRng() override { m_pixelSeed = m_globalRng(); }

  void CommitScene() override
  {
    if (m_scene)
      m_scene->Commit();
  }

  void ClearScene() override { m_scene.reset(); }

  void SetCameraPosition(float x, float y, float z) override { m_camera_position = BvhVec3(x, y, z); }

  void SetCameraRotation(float altitude, float azimuth) override
  {
    m_camera_altitude = altitude;

    m_camera_azimuth = azimuth;
  }

  void SetPerspective(float fovy, float aspect, float near, float far) override
  {
    m_perspective = Perspective{ fovy, aspect, near, far };
  }

  void SetSamplesPerPixel(int spp) override { m_samplesPerPixel = spp; }

private:
  using BvhVec3 = bvh::Vector<float, 3>;

  using BvhVec4 = bvh::Vector<float, 4>;

  using BvhRay = bvh::Ray<float>;

  template<typename Rng>
  BvhVec3 RenderPixel(int x, int y, Rng& rng) const
  {
    const float fov = std::tan(m_perspective.fovy / 2.0f);

    const float aspect = m_perspective.aspect;

    std::uniform_real_distribution<float> dist(0, 1);

    BvhVec3 totalColor(0, 0, 0);

    for (int i = 0; i < m_samplesPerPixel; i++) {

      const float xNDC = ((((x - m_viewport.x) + dist(rng)) * 2.0f) / m_viewport.w) - 1.0f;
      const float yNDC = ((((y - m_viewport.y) + dist(rng)) * 2.0f) / m_viewport.h) - 1.0f;

      const BvhVec3 rayDir = BvhVec3(xNDC * fov * aspect, yNDC * fov, -1.0f);

      const BvhRay ray(m_camera_position, rayDir);

      const int initDepth = 0;

      totalColor = totalColor + Trace(ray, rng, initDepth);
    }

    return totalColor * (1.0f / m_samplesPerPixel);
  }

  template<typename Rng>
  BvhVec3 Trace(const BvhRay& ray, Rng& rng, int depth) const
  {
    if (depth >= m_maxDepth)
      return BvhVec3(0, 0, 0);

    const auto hit = m_scene->FindClosestHit(ray);
    if (!hit)
      return OnMiss(ray);

    const BvhVec3 hitPos = ray.origin + (ray.direction * hit.distance);

    const BvhVec3 rayDir = SampleHemisphere(hit.normal, rng);

    const BvhVec3 albedo(0.8, 0.8, 0.8);

    const float shadowBias = 0.000001;

    return albedo * Trace(BvhRay(hitPos, rayDir, shadowBias), rng, depth + 1);
  }

  BvhVec3 OnMiss(const BvhRay& ray) const
  {
    const BvhVec3 up(0, 1, 0);
    const BvhVec3 lo(1.0, 1.0, 1.0);
    const BvhVec3 hi(0.5, 0.7, 1.0);
    const float level = (bvh::dot(ray.direction, up) + 1.0f) * 0.5f;
    return lo + ((hi - lo) * level);
  }

  std::shared_ptr<Geometry> DrawTriangles(int first, int count, const float* attrib, int size)
  {
    (void)first;
    (void)count;
    (void)attrib;
    (void)size;
    // TODO
    return nullptr;
  }

  std::shared_ptr<Geometry> DrawSpheres(int first, int count, const float* attrib, int size)
  {
    if (((first + count) * sizeof(float) * 4) > size_t(size)) {
      Error("Sphere attribute buffer is not large enough.");
      return nullptr;
    }

    std::shared_ptr<SphereGeometry> geom(new SphereGeometry(count));

    for (int i = 0; i < count; i++) {

      const float* vert = attrib + (i * 4);

      geom->SetCenter(i, BvhVec3(vert[0], vert[1], vert[2]));

      geom->SetRadius(i, vert[3]);
    }

    geom->Commit();

    return geom;
  }

  template<typename Rng>
  static BvhVec3 SampleHemisphere(const BvhVec3& normal, Rng& rng)
  {
    std::uniform_real_distribution<float> dist(-1, 1);

    for (int i = 0; i < 32; i++) {
      const BvhVec3 v(dist(rng), dist(rng), dist(rng));
      if ((bvh::dot(v, normal) >= 0.0f) && (bvh::dot(v, v) <= 1.0f))
        return bvh::normalize(v);
    }

    return normal;
  }

  void Error(const char* msg)
  {
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

  Viewport m_viewport;

  BvhVec3 m_camera_position = BvhVec3(0, 0, 0);

  float m_camera_altitude = 0;

  float m_camera_azimuth = 0;

  Perspective m_perspective;

  std::mt19937 m_globalRng{ 1234 };

  int m_pixelSeed = m_globalRng();

  int m_maxDepth = 3;

  int m_samplesPerPixel = 16;
};

} // namespace

std::unique_ptr<Context>
Context::Create()
{
  return std::unique_ptr<Context>(new ContextImpl);
}

} // namespace RT
