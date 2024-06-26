/**
 * @file QuadratureRule.hpp
 * @author Gianni Absillis (gabsill@ncsu.edu)
 * @brief Quadrature Rule Abstract Definition
 * @version 0.1
 * @date 2022-01-20
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once
#include <Numtool/point.hpp>

namespace iceicle {
   
    /**
     * @brief a QuadraturePoint which contains the abscisse and the quadrature weight
     * @tparam T the floating point type
     * @pram ndim the number of dimensions for the point
     */
    template<typename T, int ndim>
    struct QuadraturePoint {
        MATH::GEOMETRY::Point<T, ndim> abscisse;
        T weight;
    };

    /**
     * @brief An abstract definition of a quadrature rule
     * Quadrature rules will be defined on a reference domain
     * Quadrature rules provide quadrature points \xi_{i,g}
     * - these points are used for evaluation of the function being integrated
     * - these points are provided in the reference domain
     * 
     * Quadrature rules provide quadrature weights for integration w_g
     * 
     * an integration of the function f(\xi_i) is \sum\limits_g f(\xi_{i, g}) w_g
     * 
     * Quadrature rules also provide quadrature points and weights for the reference trace space
     * 
     * @tparam T the floating point type
     * @tparam IDX the index type
     * @tparam ndim the number of dimensions
     */
    template<typename T, typename IDX, int ndim>
    class QuadratureRule {
        public:
        /**
         * @brief Gets the number of quadrature points
         * 
         * @return int the number of quadrature points
         */
        virtual
        int npoints() const = 0;

        /**
         * @brief get the ipointh QuadraturePoint which is a struct that
         * contains the quadrature abscisse and the quadrature weight
         * @param ipoint the point index
         * @return the QuadraturePoint
         */
        virtual
        const QuadraturePoint<T, ndim> &getPoint(int ipoint) const  = 0;

        /** @brief get the quadrature point at that index (see getPoint) */
        inline
        const QuadraturePoint<T, ndim> &operator[](int ipoint) const { return getPoint(ipoint); };

        virtual 
        ~QuadratureRule() = default;
    };
}
