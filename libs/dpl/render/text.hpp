#pragma once

#include "general.hpp"

#include <glm/glm.hpp>
#include <dpl/static_vector.hpp>

namespace dpl::render
{
  class image_t
  {
    std::uint32_t w, h;
    std::vector<uint8_t> pixels;  // 1D array to store pixel data

  public:
    image_t() = default;

    image_t(std::uint32_t w, std::uint32_t h) : w{w}, h{h} {
      pixels.resize(w*h, 0);  // NOLINT(bugprone-implicit-widening-of-multiplication-result)
    }

    auto width() const {
      return w;
    }
    
    auto height() const {
      return h;
    }

    auto stride() const {
      return width();
    }

    auto size() const {
      return stride()*height();
    }
    
    auto* data() const {
      return pixels.data();
    }

    auto& operator()(std::uint32_t x, std::uint32_t y) {
      return pixels[y*w + x];
    }
  };


  
  struct pos_texel
  {
    glm::vec3 pos;
    glm::vec2 texel;
  };
  
  struct glyph_info
  {
    int advance;

    glm::ivec2 size;
    glm::ivec2 bearing;

    glm::vec2 fpos;
    glm::vec2 fsize;
  };

  inline void insertGlyph(std::vector<pos_texel>& buffer, const glyph_info& info, fb2i pen, float z0 = 0.0f)
  {
    pen->x += info.bearing.x;
    pen->y -= info.bearing.y;
    
    buffer.push_back({
      {pen->x, pen->y, z0},
      info.fpos});
    buffer.push_back({
      {pen->x, pen->y + info.size.y, z0},
      {info.fpos.x, info.fpos.y + info.fsize.y}});
    buffer.push_back({
      {pen->x + info.size.x, pen->y + info.size.y, z0},
      {info.fpos.x + info.fsize.x, info.fpos.y + info.fsize.y}});

    buffer.push_back({
      {pen->x + info.size.x, pen->y + info.size.y, z0},
      {info.fpos.x + info.fsize.x, info.fpos.y + info.fsize.y}});
    buffer.push_back({
      {pen->x + info.size.x, pen->y, z0},
      {info.fpos.x + info.fsize.x, info.fpos.y}});
    buffer.push_back({
      {pen->x, pen->y, z0},
      info.fpos});
  }

  enum class align_t : uint8_t {
    left,
    right,
    top,
    bottom,
    center
  };
  
  struct atlas_t
  {
    image_t image;

    std::map<char, glyph_info> glyphs;

    struct
    {
      std::vector<pos_texel> data;
      
      VkBuffer buffer;
      VkDeviceMemory memory;
    } recorded;

    // static constexpr int VTK_TEXT_LEFT = 0;
    // static constexpr int VTK_TEXT_CENTERED = 1;
    // static constexpr int VTK_TEXT_BOTTOM = 2;
    // static constexpr int VTK_TEXT_RIGHT = 3;
    // static constexpr int VTK_TEXT_TOP = 4;

    vector2i get_size(const std::string& str, const int extra_advance = 0) {
      int width = 0;
      int ascent = 0;
      int descent = 0;

      for (auto c : str) {
        const auto& info = glyphs[c];
        width += info.advance + extra_advance;
        ascent = std::max(ascent, info.bearing[1]);
        descent = std::max(descent, info.size[1] - info.bearing[1]);
      }

      /* last character */
      {
        const auto& last = glyphs[str.back()];
        width += last.size[0] - (last.advance + extra_advance);
      }

      return {width, ascent + descent};
    }
    
    void record(
      std::string_view str,
      fb2i pen,
      vector_n<align_t, 2> align = {align_t::left, align_t::top},
      const int extra_advance = 0) {

      int start = recorded.data.size();                           // NOLINT(cppcoreguidelines-narrowing-conversions)
      int width = 0;
      int ascent = 0;
      int descent = 0;

      /* first character */
      pen->x -= glyphs[str.front()].bearing[0];
      
      for (auto c : str) {
        const auto& info = glyphs[c];
        insertGlyph(recorded.data, info, pen, 0.5f);
        pen->x += info.advance + extra_advance;

        width += info.advance + extra_advance;        // TODO
        ascent = std::max(ascent, info.bearing[1]);   // TODO
        descent = std::max(descent, info.size[1] - info.bearing[1]);
      }

      /* last character */
      {
        const auto& last = glyphs[str.back()];
        width += last.size[0] - (last.advance + extra_advance);
      }

      // for (int i = start; i < recorded.data.size(); ++i)
      //     recorded.data[i].pos.y() += ascent;  // NOLINT
      // return;

      if (align.x == align_t::right)
        for (int i = start; i < recorded.data.size(); ++i)
          recorded.data[i].pos.x -= width;  // NOLINT
      else if (align.x == align_t::center)
        for (int i = start; i < recorded.data.size(); ++i)
          recorded.data[i].pos.x -= width/2;  // NOLINT

      
      /* Y is downward */
      if (align.y == align_t::top)
        for (int i = start; i < recorded.data.size(); ++i)
          recorded.data[i].pos.y += ascent;  // NOLINT
      else if (align.y == align_t::center)
        for (int i = start; i < recorded.data.size(); ++i)
          recorded.data[i].pos.y += ascent - (ascent + descent)/2/*height/2*/;  // NOLINT
      else
        for (int i = start; i < recorded.data.size(); ++i)
          recorded.data[i].pos.y -= descent;  // NOLINT
    }
  };
  
}
