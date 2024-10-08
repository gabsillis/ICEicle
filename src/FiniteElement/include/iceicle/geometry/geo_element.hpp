/**
 * @file geo_element.hpp
 * @author Gianni Absillis (gabsill@ncsu.edu)
 * @brief Abstract definition for Geometric Elements
 * @date 2023-06-26
 */

#pragma once
#include "iceicle/build_config.hpp"
#include <iceicle/fe_definitions.hpp>
#include <Numtool/matrixT.hpp>
#include <Numtool/fixed_size_tensor.hpp>
#include <mdspan/mdspan.hpp>
#include <memory>

namespace iceicle {

    // the maximum dynamic element order that is generated
    static constexpr int MAX_DYNAMIC_ORDER = build_config::FESPACE_BUILD_PN;

    /// @brief Collection of attributes that will uniquely identify each Element Transformation
    template<int ndim>
    struct ElementGeometryAttributes {
        DOMAIN_TYPE domain_type;
        int order;
    };

    /// @brief Represents a transformation T : s -> x which takes reference space to physical space
    /// As a function table 
    ///
    /// References:
    ///  (Guermond FE1) : Ern, Alexandre and Guermond, Jean-Luc Finite Elements I: Approximation and Interpolation
    ///  Texts in Applied Mathematics, Springer
    ///
    /// @tparam T the real number type 
    /// @tparam IDX the index type
    /// @tparam ndim the number of dimensions
    template<class T, class IDX, int ndim>
    struct ElementTransformation {

        // type aliases
        using Point = MATH::GEOMETRY::Point<T, ndim>;
        using HessianType = NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim, ndim>;
        using JacobianType = NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim>;

        // =======================
        // = Element Information =
        // =======================

        /// @brief the domain type
        DOMAIN_TYPE domain_type;

        /// @brief the order of the polynomial approximation of the transformation function
        int order;

        /// @brief the number of nodes
        int nnode;

        /// @brief the number of faces 
        int nfac;

        // ===================
        // = Node Operations =
        // ===================
        //
        /// @brief get the element coordinates from the global node coordinates and the node indices 
        ///
        /// @param [in] coord the global node coordinates array on this process
        /// @param [in] nodes the node indices 
        /// @return vector of coordinates for the nodes of this element
        std::vector<Point> (*get_el_coord)(const NodeArray<T, ndim>& coord, const IDX* nodes) = nullptr;

        // =============================
        // = Coordinate Transformation =
        // =============================

        /// @brief transform a reference domain point to the physical domain
        ///
        /// @param [in] el_coord array of the coordinates of each node, in order, for the element 
        /// @param [in] pt_ref the reference domain 
        /// @return the Point in the physical domain
        Point (*transform)(std::span<Point> el_coord, const Point& pt_ref) = nullptr;

        /// @brief get the Jacobian matrix of the transformation
        /// J = \frac{\partial T(s)}{\partial s} = \frac{\partial x}[\partial \xi}
        ///
        /// @param [in] el_coord array of the coordinates of each node, in order, for the element 
        /// @param [in] xi the position in the reference domain at which to calculate the Jacobian
        /// @return the Jacobian matrix
        JacobianType (*jacobian)(std::span<Point> el_coord, const Point& xi) = nullptr;

        /// @brief get the Hessian of the transformation
        /// H_{kij} = \frac{\partial T(s)_k}{\partial s_i \partial s_j} 
        ///         = \frac{\partial x_k}{\partial \xi_i \partial \xi_j}
        ///
        /// @param [in] el_coord array of the coordinates of each node, in order, for the element 
        /// @param [in] xi the position in the reference domain at which to calculate the hessian
        /// @return the Hessian in tensor form indexed [k][i][j] as described above
        HessianType (*hessian)(std::span<Point> el_coord, const Point& xi) = nullptr;

        // ====================================
        // = Face Transformation Connectivity =
        // ====================================

        /// @brief get the domain type of the face at the given face number
        DOMAIN_TYPE (*face_domain_type)(int face_number) = nullptr;

        /// @brief get the number of vertices in a face
        int (*n_face_vert)(int face_number) = nullptr;

        /// @brief get the vertex indices on the face
        /// NOTE: These vertices must be in the same order as if get_element_vert() 
        /// was called on the transformation corresponding to the face
        std::vector<IDX> (*get_face_vert)(int face_number, std::span<IDX> el_nodes) = nullptr;

        /// @brief get the number of nodes on the face 
        /// @param face_number the face number
        int (*n_face_nodes)(int face_number) = nullptr;

        /// @brief get the node indices on the face
        ///
        /// NOTE: Nodes are all the points defining geometry (vertices are endpoints)
        ///
        /// NOTE: These nodes must be in the same order as if get_nodes
        /// was called on the transformation corresponding to the face
        std::vector<IDX> (*get_face_nodes)(int face_number, std::span<IDX> el_nodes) = nullptr;

        /// @brief get the face number of the given vertices 
        /// @return the face number of the face with the given vertices or -1 if not found
        int (*get_face_nr) (std::span<IDX> vert_fac, std::span<IDX> el_nodes) = nullptr;


        /**
         * @brief calculate the centroid in the reference domain 
         * @return the centroid in the reference domain 
         */
        inline constexpr
        MATH::GEOMETRY::Point<T, ndim> centroid_ref() const {

            MATH::GEOMETRY::Point<T, ndim> refpt;
            // get the centroid in the reference domain based on domain type
            switch(domain_type){
                case DOMAIN_TYPE::HYPERCUBE:
                    for(int idim = 0; idim < ndim; ++idim) 
                        { refpt[idim] = 0.0; }
                    break;

                case DOMAIN_TYPE::SIMPLEX:
                    for(int idim = 0; idim < ndim; ++idim) 
                        { refpt[idim] = 1.0 / 3.0; }
                    break;

                default: // WARNING: use 0 as centroid for default
                    for(int idim = 0; idim < ndim; ++idim) 
                        { refpt[idim] = 0.0; }
                    break;
            }
            return refpt;
        }

        /**
         * @brief calculate the centroid in the physical domain 
         * @param [in] el_coord array of the coordinates of each node, in order, for the element 
         */
        inline constexpr
        MATH::GEOMETRY::Point<T, ndim> centroid( std::span<Point> el_coord ) const {
            return transform(el_coord, centroid_ref());
        }

    };
    
    /**
     * @brief A Geometric element
     * Contains information and methods for geometric description of elements
     *
     * @tparam T the floating point type
     * @tparam IDX the index type for large lists
     * @tparam ndim the number of dimensions
     */
    template<typename T, typename IDX, int ndim> 
    class GeometricElement {
        public:

        // type aliases
        using Point = MATH::GEOMETRY::Point<T, ndim>;
        using HessianType = NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim, ndim>;
        using JacobianType = NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim>;
        using value_type = T;
        using index_type = IDX;


        // =============================
        // = Coordinate Transformation =
        // =============================

        /**
         * @brief transform from the reference domain to the physical domain
         * @param [in] node_coords the coordinates of all the nodes
         * @param [in] pt_ref the point in the refernce domain
         * @param [out] pt_phys the point in the physical domain
         */
        virtual 
        void transform(
            NodeArray<T, ndim> &node_coords,
            const Point &pt_ref,
            Point &pt_phys
        ) const = 0;

        /**
         * @brief get the Jacobian matrix of the transformation
         * J = \frac{\partial T(s)}{\partial s} = \frac{\partial x}[\partial \xi}
         * @param [in] node_coords the coordinates of all the nodes
         * @param [in] xi the position in the reference domain at which to calculate the Jacobian
         * @return the Jacobian matrix
         */
        virtual
        auto Jacobian(
            NodeArray<T, ndim> &node_coords,
            const Point &xi
        ) const -> JacobianType = 0;

        /**
         * @brief get the Hessian of the transformation
         * H_{kij} = \frac{\partial T(s)_k}{\partial s_i \partial s_j} 
         *         = \frac{\partial x_k}{\partial \xi_i \partial \xi_j}
         * @param [in] node_coords the coordinates of all the nodes
         * @param [in] xi the position in the reference domain at which to calculate the hessian
         * @return the Hessian in tensor form indexed [k][i][j] as described above
         */
        virtual
        auto Hessian(
            NodeArray<T, ndim> &node_coords,
            const Point &xi
        ) const -> HessianType = 0;

        // ===============
        // = Node Access =
        // ===============

        /**
         * @brief Get the number of nodes
         * 
         * @return int the number of nodes for this element
         */
        virtual
        constexpr int n_nodes() const = 0;

        /**
         * @brief Get the nodes array for this element
         * 
         * @return pointer to the start of the array of nodes
         */
        virtual
        const IDX *nodes() const = 0;


        /**
         * @brief get the array of node indices as a span 
         *
         * This array should be garuanteed to be in the same order as the reference degres of freedom
         * for the reference domain corresponding to this face
         * This is so nodal basis functions can be mapped to the global node dofs
         *
         * @return the node indices in the mesh nodes array
         */
        std::span<const IDX> nodes_span() const {
            return std::span{nodes(), static_cast<std::size_t>(n_nodes())};
        }

        // =====================
        // = Domain Definition =
        // =====================

        /**
         * @brief get the domain type 
         * @return what refernce domain this maps to 
         */
        virtual
        constexpr DOMAIN_TYPE domain_type() const noexcept = 0;

        /**
         * @brief get the polynomial order of the geometry definition 
         * this is used to map to output
         */
        virtual 
        constexpr int geometry_order() const noexcept = 0;

        // =====================
        // = Face Connectivity =
        // =====================
        
        /// @brief get the number of faces
        virtual
        auto n_faces() const -> int = 0;

        /// @brief get the domain type of the face at the given face number
        virtual 
        auto face_domain_type(int face_number) const -> DOMAIN_TYPE = 0;

        /// @brief get the number of vertices in a face
        virtual 
        auto n_face_vert(
            int face_number /// [in] the face number
        ) const -> int = 0;

        /// @brief get the vertex indices on the face
        /// NOTE: These vertices must be in the same order as if get_element_vert() 
        /// was called on the transformation corresponding to the face
        virtual 
        auto get_face_vert(
            int face_number,      /// [in] the face number
            index_type* vert_fac  /// [out] the indices of the vertices of the given face
        ) const -> void = 0;

        /// @brief get the number of nodes on the face 
        /// @param face_number the face number
        virtual
        auto n_face_nodes(int face_number) const -> int = 0;

        /// @brief get the node indices on the face
        ///
        /// NOTE: Nodes are all the points defining geometry (vertices are endpoints)
        ///
        /// NOTE: These vertices must be in the same order as if get_nodes
        /// was called on the transformation corresponding to the face
        virtual 
        auto get_face_nodes(
            int face_number,      /// [in] the face number
            index_type* nodes_fac /// [out] the indices of the nodes of the given face
        ) const -> void = 0;

        /// @brief get the face number of the given vertices 
        /// @return the face number of the face with the given vertices
        virtual 
        auto get_face_nr(
            index_type* vert_fac /// [in] the indices of the vertices of the given face
        ) const -> int = 0;

        // =============
        // = Geometric =
        // =============

        /**
         * @brief calculate the centroid in the reference domain 
         * @return the centroid in the reference domain 
         */
        virtual 
        MATH::GEOMETRY::Point<T, ndim> centroid_ref() const {

            MATH::GEOMETRY::Point<T, ndim> refpt;
            // get the centroid in the reference domain based on domain type
            switch(domain_type()){
                case DOMAIN_TYPE::HYPERCUBE:
                    for(int idim = 0; idim < ndim; ++idim) 
                        { refpt[idim] = 0.0; }
                    break;

                case DOMAIN_TYPE::SIMPLEX:
                    for(int idim = 0; idim < ndim; ++idim) 
                        { refpt[idim] = 1.0 / 3.0; }
                    break;

                default: // WARNING: use 0 as centroid for default
                    for(int idim = 0; idim < ndim; ++idim) 
                        { refpt[idim] = 0.0; }
                    break;
            }
            return refpt;
        }

        /**
         * @brief calculate the centroid in the physical domain 
         * @param [in] coord the node coordinates arrray
         * @return the centroid in the physical domain
         */
        virtual 
        MATH::GEOMETRY::Point<T, ndim> centroid(
            NodeArray<T, ndim> &node_coords
        ) const {
            MATH::GEOMETRY::Point<T, ndim> physpt;
            
            transform(node_coords, centroid_ref(), physpt);
            return physpt;
        }

        /**
         * @brief given surface nodes 
         * find interior nodes according to their barycentric weights
         */
        virtual 
        auto regularize_interior_nodes(
            NodeArray<T, ndim>& coord /// [in/out] the node coordinates array 
        ) const -> void = 0;

        // ===========
        // = Utility =
        // ===========

        /// @brief clone this geometric element (create a copy)
        virtual 
        auto clone() const -> std::unique_ptr<GeometricElement<T, IDX, ndim>> = 0;

        /**
         * @brief virtual destructor
         */
        virtual
        ~GeometricElement() = default;
    };
}
