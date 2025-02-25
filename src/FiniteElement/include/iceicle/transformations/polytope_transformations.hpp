#pragma once
#include "Numtool/fixed_size_tensor.hpp"
#include "Numtool/point.hpp"
#include "iceicle/linalg/linalg_utils.hpp"
#include "iceicle/basis/tensor_product.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iceicle/bitset.hpp>

namespace iceicle {
 
  /// @brief namespace for polytope implementation inspired by FEMPAR
  /// Badia et al 2018
  namespace polytope {
    
    /// @brief code for a topology
    /// NOTE: indexing operator[] goes from right to left so this will match codes from FEMPAR
    ///
    /// @tparam ndim the number of dimensions
    template<std::size_t ndim>
    using tcode = bitset<ndim>;

    /// @brief code for an extrusion 
    /// An extrusion represents a marker of if the given dimension of the topology is extruded
    /// in the geometry being considered 
    /// an extrusion of all ones is the full domain of the topology 
    /// all other extrusions refer to facets
    /// @tparam ndim the number of dimensions
    template<std::size_t ndim>
    using ecode = bitset<ndim>;

    template<std::size_t ndim>
    static constexpr ecode<ndim> full_extrusion = ~ecode<ndim>{0};

    /// @brief code for a vertex 
    /// @tparam ndim the number of dimensions
    template<std::size_t ndim>
    using vcode = bitset<ndim>;

    /// @brief code for a prismatic extrusion
    static constexpr bool prism_ext = 1;

    /// @brief code for a simplical extrusion 
    static constexpr bool simpl_ext = 0;

    // ======================
    // = geo_code utilities =
    // ======================

    /// @brief convert a vertex code to a point in space 
    /// @tparam T the floating point type 
    /// @tparam ndim the number of dimensions
    template<class T, std::size_t ndim>
    constexpr 
    auto vcode_to_point(vcode<ndim> v) noexcept -> MATH::GEOMETRY::Point<T, (int) ndim>
    {
      MATH::GEOMETRY::Point<T, (int) ndim> pt{};
      for(int idim = 0; idim < ndim; ++idim){
        if(v[idim] == 0){
          pt[idim] = static_cast<T>(0.0);
        } else {
          pt[idim] = static_cast<T>(1.0);
        }
      }
    }

    /// @brief get the number of dimensions of given code
    template<std::size_t ndim>
    constexpr
    auto get_ndim(bitset<ndim> code) -> int {
      return ndim;
    }

    template<class T>
    concept geo_code = requires(T code) {
      get_ndim(code);
    };

    template<geo_code auto c>
    struct static_geo_code{};

    /// @brief get the tcode of the polytope that forms the "base" of an extruded polytope
    /// i.e the bottom square of a pyramid
    template< std::size_t ndim >
    [[nodiscard]] inline constexpr 
    auto base_tcode(bitset<ndim> tcode)
    { 
      if constexpr (ndim == 0) return bitset<0>{};
      else return bitset<ndim - 1>(tcode.to_ullong()); 
    }

    /// @brief get the tcode of the polytope 
    /// that results from a prismtic extrusion of the argument
    /// @param t the polytope to extrude 
    /// @return a prismatic extrusion of t
    template< std::size_t ndim >
    [[nodiscard]] inline constexpr 
    auto prismatic_extrude(tcode<ndim> t) 
    -> tcode<ndim + 1> 
    {
      tcode<ndim + 1> extruded(t.to_ullong());
      extruded[ndim] = prism_ext;
      return extruded;
    }

    /// @brief get the tcode of the polytope 
    /// that results from a simplical extrusion of the argument
    /// @param t the polytope to extrude 
    /// @return a simplical extrusion of t
    template< std::size_t ndim >
    [[nodiscard]] inline constexpr 
    auto simplical_extrude(tcode<ndim> t)
    -> tcode<ndim + 1> 
    {
      tcode<ndim + 1> extruded(t.to_ullong());
      extruded[ndim] = simpl_ext;
      return extruded;
    }

    // ============
    // = Vertices =
    // ============

    /// @brief get the number of vertices for the given polytope 
    /// @param t the code that defines the polytope domain 
    /// @return the number of vertices
    template<std::size_t ndim>
    constexpr
    auto n_vert(tcode<ndim> t) noexcept -> std::size_t 
    {
      // special case (we will represent 0 dimensional objects with a single vertex)
      if constexpr(ndim == 0) { return 1; }
      else {
        std::size_t nvert = 2;
        for(auto idim = 1; idim < ndim; ++idim){
          //prism extrusions double the number of vertices 
          //while simplical extrusions converge to a single vertex at alpha = 1
          if(t[idim] == prism_ext) nvert *= 2;
          else nvert += 1;
        }
        return nvert;
      }
    }

    /// @brief get the indices of vertices corresponding to the orientation axis 
    /// This is (in order)
    /// - The vertex at the origin 
    /// - The vertex one unit in the x-direction from the origin 
    /// - The vertex one unit in the y-direction from the origin 
    /// - ... 
    template<std::size_t ndim>
    [[nodiscard]] inline constexpr 
    auto orient_axis_vert(tcode<ndim> t) noexcept 
    -> std::array<std::size_t, ndim + 1>
    {
      std::array<std::size_t, ndim + 1> orient_axis;
      orient_axis[0] = 0;
      orient_axis[1] = 1;
      for(std::size_t idim = 1; idim < ndim; ++idim){
        if(t[idim] == simpl_ext){
          orient_axis[idim + 1] = orient_axis[idim] + 1;
        } else {
          orient_axis[idim + 1] = orient_axis[idim] * 2;

        }
      }
      return orient_axis;
    }

    template<geo_code auto t>
    using vertex_list = std::array<vcode<get_ndim(t)>, n_vert(t)>;

    /// @brief generate the list of vertices for a given topology 
    /// @tparam the bitcode for the topology
    /// @return a list of bitcodes for each vertex 
    /// The bitcodes correspond to the given coordinate being either 0.0 or 1.0;
    template<geo_code auto t>
    constexpr
    auto gen_vert() noexcept -> vertex_list<t>
    {
      static constexpr int ndim = get_ndim(t);
      vertex_list<t> vertices;
      std::ranges::fill(vertices, vcode<ndim>{0});

      std::size_t nvert_current = 1;
      for(int idim = 0; idim < ndim; ++idim){
        if(t[idim] == simpl_ext){
          // simplex extrusion domain comes to a single point
          vertices[nvert_current] = vcode<ndim>{static_cast<unsigned long long>(std::pow(2, idim))};
          ++nvert_current;
        } else {
          // prismatic extrusion extrudes all the vertices of the current domain
          std::copy_n(vertices.begin(), nvert_current, vertices.begin() + nvert_current);
          for(int ivert = nvert_current; ivert < 2 * nvert_current; ++ivert){
            vertices[ivert][idim] = 1;
          }
          nvert_current *= 2;
        }
      }
      return vertices;
    }

    // ==========
    // = Facets =
    // ==========

    /// @brief get the number of facets of a given dimension of a topology
    /// @param t the topology code 
    /// @param idim the dimension of the facets 
    ///       (0 is vertices, ndim is the volume)
    template<std::size_t idim>
    [[nodiscard]] inline constexpr
    auto n_facets(geo_code auto t) 
    -> std::size_t 
    {
      if constexpr (idim == 0) return n_vert(t);
      else if (idim == get_ndim(t)) return 1;
      else if (idim > get_ndim(t)) return 0;
      else {
        auto t_base{base_tcode(t)};

        std::size_t facet_count = n_facets<idim>(t_base);
        // count the number of facets on the base
        if(t[get_ndim(t) - 1] == prism_ext)
          facet_count *= 2;

        // count the number of facets in the extrusion
        // this is the number of idim - 1 facets of the base shape
        facet_count += n_facets<idim - 1>(t_base);
        return facet_count;
      }
    }


    template<std::size_t idim>
    struct facet {
      tcode<idim> t;
      std::vector<std::size_t> vertex_indices;
    };

    /// @brief generate the topology and vertex indices of all of the 
    /// idim-dimensional facets of the topology t 
    /// i.e idim = 1 on a hex toploogy will give all the lines
    /// that comprise that hex
    ///
    /// The vertex indices of the facets should be ordered such that 
    /// the orientation axis of the faces ( ndim-1 facets of a topology )
    /// generate normal vectors that are outward with respect to the topology
    template< std::size_t idim >
    [[nodiscard]] inline constexpr 
    auto get_facets(geo_code auto t)
    -> std::vector< facet<idim> >
    {
      if constexpr (idim > decltype(t)::static_extent()){
        return std::vector< facet<idim> >{};
      }
      if constexpr (idim == decltype(t)::static_extent()){
        // the desired facet is the entire polytope
        std::vector<std::size_t> vertex_indices(n_vert(t));
        std::iota(vertex_indices.begin(), vertex_indices.end(), 0);
        return std::vector{facet{t, vertex_indices}};
      } else if constexpr(idim <= 0){
        // the desired facets are all the vertices
        std::vector< facet<idim> > facets{};
        for(std::size_t idx = 0; idx < n_vert(t); ++idx){
          facets.push_back(facet{tcode<0>{}, std::vector<std::size_t>{idx}});
        }
        return facets;
      } else {
        std::vector< facet<idim> > facets{};

        // first get the vertices of the facets on the base 
        auto base_facets = get_facets<idim>(base_tcode(t));
        facets.insert(facets.end(), base_facets.begin(), base_facets.end());

        // then, if prismatic extrusion, get the facets on the mirror of the base
        // these indices are offsest by the number of vertices in the base
        if(t[get_ndim(t) - 1] == prism_ext){
          auto top_facets = get_facets<idim>(base_tcode(t));
          for(facet<idim>& f : top_facets){
            for(std::size_t& idx : f.vertex_indices){
              idx += n_vert(base_tcode(t));
            }
          }
          facets.insert(facets.end(), top_facets.begin(), top_facets.end());
        }

        // then get the facets generated by the extrusion
        auto base_subfacets = get_facets<idim - 1>(base_tcode(t));
        if(t[get_ndim(t) - 1] == prism_ext){
          for(const facet<idim - 1>& f : base_subfacets){
            std::vector<std::size_t> new_vertices = f.vertex_indices;
            for(std::size_t v : f.vertex_indices){
              new_vertices.push_back(n_vert(base_tcode(t)) + v);
            }
            tcode<idim> extruded{prismatic_extrude(f.t)};
            facets.push_back(facet<idim>{extruded, new_vertices});
          }
        } else {
          for(const facet<idim - 1>& f : base_subfacets){
            std::vector<std::size_t> new_vertices = f.vertex_indices;
            new_vertices.push_back(n_vert(base_tcode(t)));
            tcode<idim> extruded{simplical_extrude(f.t)};
            facets.push_back(facet<idim>{extruded, new_vertices});
          }
        }

        return facets;
      }
    }

    /// @brief an extrusion is said to have even parity if the sign for the 
    /// wedge product of basis vectors that form the orientation definition of the extrusion 
    ///
    /// e.g (e, v) = (110, 110)
    /// the orientation is defined by jhat /\ khat (y and z bits are set to 1 in e)
    /// the y coordinate of v == 1, therefore the j direciton of extrusion is negative (to be interior)
    /// the z coordinate of v == 1, therefore the k direction of extrusion is negative 
    /// -jhat /\ -khat = + (jhat /\ khat) therefore the parity is true
    ///
    /// is positive (returns true)
    /// if this is negative returns false
    template<std::size_t ndim>
    constexpr 
    auto extrusion_parity(ecode<ndim> e, vcode<ndim> v) -> bool 
    { return ( (e.count() - (e & ~v).count()) % 2 ) == 0; }

    /// @brief get the parity as defined above in extrusion_parity() 
    /// for the hodge dual of a given extrusion from a vertex
    template<std::size_t ndim>
    constexpr
    auto hodge_extrusion_parity(ecode<ndim> e, vcode<ndim> v) -> bool 
    {
      //TODO: probably faster way to directly compute levi civita based on e

      // set up the levi civita tensor for the hodge dual of the extrusion
      std::array<std::size_t, ndim> lc_indices{};
      int iindex = 0;

      // all dimensions in order where e == 1
      for(int idim = 0; idim < ndim; ++idim){
        if(e[idim] == 1) {
          lc_indices[iindex] = idim;
          iindex++;
        }
      }

      // all dimension in order where e == 0
      for(int idim = 0; idim < ndim; ++idim){
        if(e[idim] == 0) {
          lc_indices[iindex] = idim;
          iindex++;
        }
      }

      return extrusion_parity(e, v) ^
        (NUMTOOL::TENSOR::FIXED_SIZE::levi_civita<int, ndim>.list_index(lc_indices.data()) == 1);
    }

    /// @brief given a topology and extrusion code -- get the number of vertices 
    /// @param t the topology code 
    /// @param e the extrusion code
    template<std::size_t ndim>
    [[nodiscard]] inline constexpr 
    auto n_vert(tcode<ndim> t, ecode<ndim> e)
    -> std::size_t {
      std::size_t nvert = 1;
      for(int idim = 0; idim < ndim; ++idim){
        if(e[idim] != 0){
          if(t[idim] == prism_ext) nvert *= 2;
          else nvert += 1;
        }
      }
      return nvert;
    }

    // @brief get the topology of a given extrusion e of a topology t 
    // @tparam t the toplogy 
    // @tparam e the extrusion 
    // @return the toplogy of the e extrusion of t
    template<geo_code auto t, geo_code auto e>
    [[nodiscard]] inline constexpr
    auto extrusion_topology()
    -> tcode<e.count()> 
    {
      std::size_t jdim = 0;
      tcode<e.count()> ext_t{};
      for(std::size_t idim = 0; idim < get_ndim(t); ++idim)
        if(e[idim] == 0) ext_t[jdim++] = t[idim];
      return ext_t;
    }

    [[nodiscard]] inline constexpr
    auto validate_facet(geo_code auto t, geo_code auto e, geo_code auto v)
    -> bool {
      // prevent ambiguous facet definitions
      int ndim = get_ndim(t);
      for(int idim = 0; idim < ndim; ++idim){
        
      }
      return true;
    }

    template<geo_code auto t, geo_code auto e, geo_code auto v>
    [[nodiscard]]
    auto facet_vertices()
    -> std::array<vcode<get_ndim(t)>, n_vert(t, e)>
    {
      static constexpr int ndim = get_ndim(t);
      std::array<vcode<get_ndim(t)>, n_vert(t, e)> vertices;
      vertices[0] = v;
      std::size_t inode = 1;
      for(int idim = 0; idim < ndim; ++idim){
        if (e[idim] != 0){
          if (t[idim] == simpl_ext){
            vertices[inode] = vertices[inode - 1];
            for(int jdim = 0; jdim < ndim; ++jdim)
              vertices[inode][jdim] = 0;
            vertices[inode].flip(idim);
            ++inode;
          } else {
            std::copy_n(vertices.begin(), inode, &vertices[inode]);
            std::size_t end = 2 * inode;
            for(; inode < end; ++inode ){
              vertices[inode].flip(idim);
            }
          }
        }
      }
      return vertices;
    }

  }
}
