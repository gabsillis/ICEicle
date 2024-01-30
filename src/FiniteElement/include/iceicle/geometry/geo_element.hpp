/**
 * @file geo_element.hpp
 * @author Gianni Absillis (gabsill@ncsu.edu)
 * @brief Abstract definition for Geometric Elements
 * @date 2023-06-26
 */

#pragma once
#include "iceicle/fe_enums.hpp"
#include <iceicle/fe_function/nodal_fe_function.hpp>
#include <Numtool/point.hpp>
#include <Numtool/matrixT.hpp>
#include <Numtool/fixed_size_tensor.hpp>
#include <mdspan/mdspan.hpp>

namespace ELEMENT{

    // the maximum dynamic element order that is generated
    static constexpr int MAX_DYNAMIC_ORDER = 5;
    
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
        private:

        // namespace aliases
        using Point = MATH::GEOMETRY::Point<T, ndim>;
        using HessianType = NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim, ndim>;

        public:
        /**
         * @brief Get the number of nodes
         * 
         * @return int the number of nodes for this element
         */
        virtual
        constexpr int n_nodes() const = 0;

        /**
         * @brief get the domain type 
         * @return what refernce domain this maps to 
         */
        virtual
        constexpr FE::DOMAIN_TYPE domain_type() const noexcept = 0;

        /**
         * @brief calculate the centroid in the physical domain 
         * @param [in] coord the node coordinates arrray
         * @return the centroid in the physical domain
         */
        virtual 
        MATH::GEOMETRY::Point<T, ndim> centroid(
            FE::NodalFEFunction<T, ndim> &node_coords
        ) const {
            MATH::GEOMETRY::Point<T, ndim> refpt;
            MATH::GEOMETRY::Point<T, ndim> physpt;

            // get the centroid in the reference domain based on domain type
            switch(domain_type()){
                case FE::DOMAIN_TYPE::HYPERCUBE:
                    for(int idim = 0; idim < ndim; ++idim) 
                        { refpt[idim] = 0.0; }
                    break;

                case FE::DOMAIN_TYPE::SIMPLEX:
                    for(int idim = 0; idim < ndim; ++idim) 
                        { refpt[idim] = 1.0 / 3.0; }
                    break;

                default: // WARNING: use 0 as centroid for default
                    for(int idim = 0; idim < ndim; ++idim) 
                        { refpt[idim] = 0.0; }
                    break;
            }
            
            transform(node_coords, refpt, physpt);
            return physpt;
        }

        /**
         * @brief get the polynomial order of the geometry definition 
         * this is used to map to output
         */
        virtual 
        constexpr int geometry_order() const noexcept = 0;

        /**
         * @brief Get the nodes array for this element
         * 
         * @return std::ptrdiff_t* the array of nodes
         */
        virtual
        const IDX *nodes() const = 0;

        /**
         * @brief transform from the reference domain to the physical domain
         * @param [in] node_coords the coordinates of all the nodes
         * @param [in] pt_ref the point in the refernce domain
         * @param [out] pt_phys the point in the physical domain
         */
        virtual 
        void transform(
            FE::NodalFEFunction<T, ndim> &node_coords,
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
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim> Jacobian(
            FE::NodalFEFunction<T, ndim> &node_coords,
            const Point &xi
        ) const = 0;

        /**
         * @brief get the Hessian of the transformation
         * H_{kij} = \frac{\partial T(s)_k}{\partial s_i \partial s_j} 
         *         = \frac{\partial x_k}{\partial \xi_i \partial \xi_j}
         * @param [in] node_coords the coordinates of all the nodes
         * @param [in] xi the position in the reference domain at which to calculate the hessian
         * @return the Hessian in tensor form indexed [k][i][j] as described above
         */
        virtual
        HessianType Hessian(
            FE::NodalFEFunction<T, ndim> &node_coords,
            const Point &xi
        ) const = 0;

        /**
         * @brief virtual destructor
         */
        virtual
        ~GeometricElement() = default;
    };
}
