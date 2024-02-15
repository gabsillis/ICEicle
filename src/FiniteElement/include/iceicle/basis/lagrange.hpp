/**
 * @file lagrange.hpp
 * @author Gianni Absillis (gabsill@ncsu.edu)
 * @brief Lagrange Basis functions
 * @version 0.1
 * @date 2023-06-26
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#pragma once
#include "iceicle/fe_enums.hpp"
#include "iceicle/transformations/HypercubeElementTransformation.hpp"
#include <algorithm>
#include <iceicle/basis/basis.hpp>
#include <iceicle/transformations/SimplexElementTransformation.hpp>
#include <mdspan/mdspan.hpp>
namespace BASIS {
    
    /**
     * @brief Lagrange Basis functions on simplex elements
     * @tparam T the floating point type
     * @tparam IDX the index type
     * @tparam ndim the number of dimensions
     * @tparam Pn the polynomial order
     */
    template<typename T, typename IDX, int ndim, int Pn>
    class SimplexLagrangeBasis final : public Basis<T, ndim>{
        private:

        static inline ELEMENT::TRANSFORMATIONS::SimplexElementTransformation<T, IDX, ndim, Pn> transform;

        public:

        // ==============
        // = Basis Impl =
        // ==============

        int nbasis() const override { return transform.nnodes(); }

        constexpr FE::DOMAIN_TYPE domain_type() const noexcept override { return FE::DOMAIN_TYPE::SIMPLEX; }

        void evalBasis(const T *xi, T *Bi) const override {
            for(int inode = 0; inode < transform.nnodes(); ++inode){
                Bi[inode] = transform.shp(xi, inode);
            }
        }

        void evalGradBasis(const T *xi, T *dBidxj) const override 
        {
            auto dB = std::experimental::mdspan(dBidxj, transform.nnodes(), ndim);
            for(int inode = 0; inode < transform.nnodes(); ++inode){
                for(int jderiv = 0; jderiv < ndim; ++jderiv){
                   dB[inode, jderiv] = transform.dshp(xi, inode, jderiv);
                }
            }
        }

        bool isOrthonormal() const override { return false; }

        bool isNodal() const override { return true; }

        inline int getPolynomialOrder() const override { return Pn; }
    };


    template<typename T, typename IDX, int ndim, int Pn>
    class HypercubeLagrangeBasis final: public Basis<T, ndim> {
        static inline ELEMENT::TRANSFORMATIONS::HypercubeElementTransformation<T, IDX, ndim, Pn> transform{};

        using Point = MATH::GEOMETRY::Point<T, ndim>;
        public:

        // ==============
        // = Basis Impl =
        // ==============

        int nbasis() const override { return transform.n_nodes(); }

        constexpr FE::DOMAIN_TYPE domain_type() const noexcept override { return FE::DOMAIN_TYPE::HYPERCUBE; }

        void evalBasis(const T*xi, T *Bi) const override {
            Point xipt{};
            std::copy_n(xi, ndim, xipt.data());
            transform.fill_shp(xipt, Bi);
        }

        void evalGradBasis(const T *xi, T *dBidxj) const override {
            Point xipt{};
            std::copy_n(xi, ndim, xipt.data());
            NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, transform.n_nodes(), ndim> dBi; 
            transform.fill_deriv(xipt, dBi);
            // TODO: get rid of memmove called here
            std::copy_n(dBi.ptr(), ndim * transform.n_nodes(), dBidxj);
        }
        
        void evalHessBasis(const T *xi, T *HessianData) const override {
            transform.fill_hess(xi, HessianData);
        }

        bool isOrthonormal() const override { return false; }

        bool isNodal() const override { return true; }

        inline int getPolynomialOrder() const override { return Pn; }
    };
}
