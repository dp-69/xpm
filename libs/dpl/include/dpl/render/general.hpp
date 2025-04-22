#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <tuple>

namespace dpl::render
{
  struct ndc2f
  {
    glm::vec2 v;

    auto& operator*() const { return v; }
    
    auto* operator->() { return &v; }
    auto* operator->() const { return &v; }

    ndc2f() = default;

    ndc2f(float x, float y) : v{x, y} {}

    ndc2f(const glm::mat4& m, const glm::vec3& pos) {  // NOLINT(cppcoreguidelines-pro-type-member-init)
      auto p = m*glm::vec4{pos.x, pos.y, pos.z, 1.0f};
      v.x = p.x/p.w;
      v.y = p.y/p.w;
    }
  };

  struct ndc3f
  {
    glm::vec3 v;

    auto& operator*() const { return v; }
    
    auto* operator->() { return &v; }
    auto* operator->() const { return &v; }

    ndc3f(float x, float y, float z) : v{x, y, z} {}
    
    ndc3f(const glm::mat4& m, const glm::vec3& pos) {
      auto p = m*glm::vec4{pos.x, pos.y, pos.z, 1.0f};
      v.x = p.x/p.w;
      v.y = p.y/p.w;
      v.z = p.z/p.w;
    }

    ndc2f add2(const glm::vec2& rhs) const {
      return ndc2f{v.x + rhs.x, v.y + rhs.y};
    }

    ndc3f add3(const glm::vec2& rhs) const {
      return ndc3f{v.x + rhs.x, v.y + rhs.y, 0.0f};
    }
  };


  template<typename T>
  struct fb2t
  {
    glm::vec<2, T> v;

    auto& operator*() const { return v; }
    
    auto* operator->() { return &v; }

    void set(const ndc2f& ndc, const glm::vec2& half_extent) {
      v.x = (ndc->x + 1)*half_extent.x;
      v.y = (ndc->y + 1)*half_extent.y;
    }

    fb2t() = default;

    fb2t(T x, T y) : v{x, y} {}

    fb2t(const ndc2f& ndc, const glm::vec2& half_extent, std::true_type /* round */) {
      v.x = std::round((ndc->x + 1)*half_extent.x);
      v.y = std::round((ndc->y + 1)*half_extent.y);
    }
  };

  using fb2f = fb2t<float>;
  using fb2i = fb2t<int>;
  


  


  

  struct magnification_t
  {
    float factor = 0;

    void reset() {
      factor = 0;
    }
  };
  
  struct rotation_t
  {
    float angle = 0;
    glm::vec3 axis{0, 1, 0};

    auto camera_rotation() const {
      return rotate(glm::mat4{1}, angle, axis);
    }

    void reset() {
      angle = 0;
      axis = {0, 1, 0};
    }
  };

  struct translation_t
  {
    glm::vec3 vec{0, 0, 0};

    void reset() {
      vec = {0, 0, 0};
    }
  };

  struct camera_t
  {

    glm::vec3 eye   {7.94687, -2.84443, -2.33016};
    glm::vec3 center{1.46804, 0.491919, 0.117398};
    glm::vec3 up    {-0.268466, -0.851608, 0.45021};
    
    // glm::vec3 eye   {0,  0, -4};
    // glm::vec3 center{0,  0,  0};
    // glm::vec3 up    {0, -1,  0};

    
    // glm::vec3 eye   {-4,  0.5f, 0.5f}; //{0,  0, -4};
    // glm::vec3 center{0.5f}; //{0,  0,  0};
    // glm::vec3 up    {-0.122935, -0.991247, -0.0481162};// {0, -1,  0};
    
    // glm::vec3 eye   {-3.9944, -0.173819, -1.0578}; //{0,  0, -4};
    // glm::vec3 center{0.298674, 0.212657, 0.422827}; //{0,  0,  0};
    // glm::vec3 up    {-0.122935, -0.991247, -0.0481162};// {0, -1,  0};

    auto dir() const {
      return center - eye;
    }

    auto right() const {
      return cross(normalize(dir()), up);
    }

    auto rotated(const rotation_t& r) const {
      auto mat = r.camera_rotation();  // return std::make_tuple(eye, up);

      using namespace glm;
      return std::make_tuple(
        center - vec3(mat*vec4(dir(), 1)),
        vec3(mat*vec4(up, 1)));
    }

    void update_magnification(magnification_t& m) {
      eye = center - (1 - m.factor)*dir();
      m.reset();
    }

    void update_translation(translation_t& t) {
      eye    += t.vec;
      center += t.vec;
      t.reset();
    }

    void update_rotation(rotation_t& r) {
      std::tie(eye, up) = rotated(r);
      r.reset();
    }
  };  







  /**
    * https://www.youtube.com/watch?v=rvJHkYnAR3w&t=261s&ab_channel=BrendanGalea
    *
    * @param eye 
    * @param dir normal
    * @param up normal
    * @return 
    */
  inline glm::mat4 viewDir(const glm::vec3& eye, const glm::vec3& dir, const glm::vec3& up) {
    glm::vec3 r{cross(dir, up)};

    glm::mat4 mat{1};
    mat[0][0] =   r.x;   mat[1][0] =   r.y;   mat[2][0] =   r.z;   mat[3][0] = -dot(r, eye);
    mat[0][1] = -up.x;   mat[1][1] = -up.y;   mat[2][1] = -up.z;   mat[3][1] =  dot(up, eye);
    mat[0][2] = dir.x;   mat[1][2] = dir.y;   mat[2][2] = dir.z;   mat[3][2] = -dot(dir, eye);

    return mat;
  }

  inline glm::mat4 viewAt(const glm::vec3& eye, const glm::vec3& center, const glm::vec3& up) {
    return viewDir(eye, normalize(center - eye), up);
  }

  /**
   * https://www.youtube.com/watch?v=U0_ONQQ5ZNM&ab_channel=BrendanGalea
   *
   * @param near can be zero
   */
  inline glm::mat4 orthographic(float left, float right, float top, float bottom, float near = 0.0f, float far = 1.0f) {
    glm::mat4 mat{1};

    mat[0][0] = 2/(right - left);
    mat[1][1] = 2/(bottom - top);
    mat[2][2] = 1/(far - near);
    mat[3][0] = (right + left)/(left - right);
    mat[3][1] = (bottom + top)/(top - bottom);
    mat[3][2] = near*-mat[2][2];

    return mat;
  }

  inline float cross_z(const glm::vec2& a, const glm::vec2& b) {
    return a.x*b.y - a.y*b.x; 
  }

  /**
   * https://www.youtube.com/watch?v=U0_ONQQ5ZNM&ab_channel=BrendanGalea
   *
   * @param near must not be zero
   */
  inline glm::mat4 perspective(float fovy, float aspect, float near, float far) {
    float invTanHalfFovy = 1/tan(fovy/2);

    glm::mat4 mat{0};

    mat[0][0] = invTanHalfFovy/aspect;
    mat[1][1] = invTanHalfFovy;
    mat[2][2] = far/(far - near);
    mat[2][3] = 1;
    mat[3][2] = near*-mat[2][2];

    return mat;
  }
}