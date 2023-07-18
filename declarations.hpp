#pragma once

#include <dpl/static_vector.hpp>

#include <boost/pending/disjoint_sets.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <iostream>
#include <fstream>



namespace xpm
{
  using v3i = dpl::vector3i;
  using v3d = dpl::vector3d;

  namespace voxel_tag
  {
    struct phase
    {
      std::uint8_t value;

      friend bool operator==(const phase& lhs, const phase& rhs) { return lhs.value == rhs.value; }
      friend bool operator!=(const phase& lhs, const phase& rhs) { return !(lhs == rhs); }

      const auto& operator*() const {
        return value;
      }
    };

    /**
     * \brief
     *   -1  : for non-pore voxels
     *   >=0 : pore node of a pore-voxel
     */
    struct velem
    {
      std::int32_t value;

      // auto& operator*() {
      //   return value;
      // }

      const auto& operator*() const {
        return value;
      }
    };


  }

  namespace presets
  {
    static inline constexpr voxel_tag::phase pore = {0};
    static inline constexpr voxel_tag::phase solid = {1};
    static inline constexpr voxel_tag::phase microporous = {2};

    static constexpr auto darcy_to_m2 = 9.869233e-13;
  }

  namespace geometric_properties
  {
    struct equilateral_triangle_properties
    {
      static constexpr double area(double r_ins = 1) {
        return 5.19615242271*r_ins*r_ins;
      }

      // k * G, k - coefficient, G - shape factor
      static constexpr double conductance(double area = 1, double viscosity = 1) {
        return 0.0288675134595*area*area/viscosity;   // = std::sqrt(3)/60 = k*G*A^2/mu for eq tri
      }
    };
  }
  
  namespace attribs {
    def_static_key(pos)
    def_static_key(r_ins)
    def_static_key(adj)
    def_static_key(length)
    def_static_key(length0)
    def_static_key(length1)
  }


  /**
   * \brief maximum node count
   */
  using idx1d_t = int32_t;
  using idx3d_t = dpl::vector_n<idx1d_t, 3>;

  using disjoint_sets = boost::disjoint_sets_with_storage<
    boost::typed_identity_property_map<idx1d_t>,
    boost::typed_identity_property_map<idx1d_t>
  >;



  // struct idx1d_expl
  // {
  //   idx1d_t value;
  // };

  // struct idx3d_expl
  // {
  //   idx3d_t value;
  // };

  // template<typename T, int n> requires std::integral<T>
  // class map_idx_t
  // {
  //   // dpl::vector_n<Type, n> dim;
  // };

  template <typename R> requires std::integral<R>
  class map_idx3_t
  {
    R x_, xy_;

  public:
    template<typename T>
    explicit map_idx3_t(dpl::vector_n<T, 3> dim)
      : x_(dim.x()), xy_(static_cast<R>(dim.x())*dim.y()) {}

    template<typename V>
    R operator()(dpl::vector_n<V, 3> v) const {
      return static_cast<R>(v.x()) + x_*v.y() + xy_*v.z();
    }

    template<typename V>
    R operator()(V x, V y, V z) const {
      return static_cast<R>(x) + x_*y + xy_*z;
    }

    template <typename I>
    auto operator[](std::integral_constant<I, 0>) { return 1; }
     
    template<typename I>
    auto operator[](std::integral_constant<I, 1>) { return x_; }

    template <typename I>
    auto operator[](std::integral_constant<I, 2>) { return xy_; }
  };


  inline auto idx_mapper(idx3d_t dim) {
    return map_idx3_t<idx1d_t>{dim};
  }













  namespace parse
  {
    struct image_dict
    {
      std::uint8_t solid;
      std::uint8_t pore;
      std::uint8_t microporous;
    };
  }


  


  


  
}
