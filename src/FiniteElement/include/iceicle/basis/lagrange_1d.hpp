/**
 * @brief 1D Lagrange basis functions
 * @author Gianni Absillis (gabsill@ncsu.edu)
 */
#pragma once
#include "Numtool/polydefs/LagrangePoly.hpp"
#include <Numtool/fixed_size_tensor.hpp>

namespace iceicle {

    /**
     * @brief Interpolation of a uniform set of Pn + 1 points from -1 to 1 
     */
    template<typename T, int Pn>
    struct UniformLagrangeInterpolation {
        template<typename T1, std::size_t... sizes>
        using Tensor = NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T1, sizes...>;

        /// @brief value type tag
        using value_type = T;

        /// @brief the number of basis functions generated by this interpolation
        static constexpr int nbasis = Pn + 1;

        /// compile-time precompute the evenly spaced interpolation points
        static constexpr Tensor<T, Pn + 1> xi_nodes = []{
        Tensor<T, Pn + 1> ret = {};
        if constexpr (Pn == 0) {
            // finite volume should recover cell center
            // for consistency
            ret[0] = 0.0; 
        } else {
            T dx = 2.0 / Pn;
            ret[0] = -1.0;
            for(int j = 1; j < Pn + 1; ++j){
            // better for numerics than j * dx
            ret[j] = ret[j - 1] + dx;
            }
        }
        return ret;
        }();

        /// compile-time precompute the lagrange polynomial denominators
        /// the barycentric weights
        static constexpr Tensor<T, Pn + 1> wj = []{
            Tensor<T, Pn + 1> ret;
            for(int j = 0; j < Pn + 1; ++j){
                ret[j] = 1.0;
                for(int k = 0; k < Pn + 1; ++k) if(k != j) {
                ret[j] *= (xi_nodes[j] - xi_nodes[k]);
                }
                // invert (this is a denominator)
                // NOTE: Berrut, Trefethen have an optimal way to compute this 
                // but we don't because its computed at compile time
                ret[j] = 1.0 / ret[j];
            }
            return ret;
        }();

        /**
        * @brief Evaluate every interpolating polynomial at the given point 
        * @param xi the point to evaluate at 
        * @return an array of all the evaluations
        */
        Tensor<T, Pn + 1> eval_all(T xi) const {
            Tensor<T, Pn + 1> Nj{};

            // finite volume case
            if constexpr (Pn == 0){
                Nj[0] = 1.0;
            } else if constexpr (Pn == 1){
                // hard code the simple cases for efficiency and reproducibility
                // (Tiny amounts of roundoff with general case)
                Nj[0] = 0.5 * (1 - xi);
                Nj[1] = 1.0 - Nj[0];
            } else {
                // run-time precompute the product of differences
                T lskip = 1; // this is the product skipping the node closest to xi
                int k; // this will be used to determine which node to skip
                for(k = 0; k < Pn; ++k){
                    // make sure k+1 isn't the closest node to xi
                    if( xi >= (xi_nodes[k] + xi_nodes[k+1])/2 ){
                    lskip *= (xi - xi_nodes[k]);
                    } else {
                    break; // k is the closest node to xi
                    }
                }
                for(int i = k + 1; i < Pn + 1; ++i){
                    lskip *= (xi - xi_nodes[i]);
                }
                T lprod = lskip * (xi - xi_nodes[k]);

                // calculate Nj 
                int j;
                for(j = 0; j < k; ++j){
                    Nj[j] = lprod * wj[j] / (xi - xi_nodes[j]);
                }
                Nj[k] = lskip * wj[k];
                for(++j; j < Pn + 1; ++j){
                    Nj[j] = lprod * wj[j] / (xi - xi_nodes[j]);
                }
            }

            return Nj;
        }

        /**
        * @brief Get the value and derivative of 
        * every interpolating polynomial at the given point 
        * @param xi the point to evaluate at 
        * @return an array of all the evaluations
        */
        void deriv_all(
            T xi,
            Tensor<T, Pn+1> &Nj,
            Tensor<T, Pn+1> &dNj
        ) const {

            // finite volume case
            if constexpr (Pn == 0){
                Nj[0] = 1.0;
                dNj[0] = 0.0;
            } else if constexpr (Pn == 1){

                // hard code the simple cases for efficiency and reproducibility
                // (Tiny amounts of roundoff with general case)
                Nj[0] = 0.5 * (1 - xi);
                Nj[1] = 1.0 - Nj[0];
                dNj[0] = -0.5;
                dNj[1] = 0.5;
            } else {
                // run-time precompute the product of differences
                T lskip = 1; // this is the product skipping the node closest to xi
                int k; // this will be used to determine which node to skip
                for(k = 0; k < Pn; ++k){
                    // make sure k+1 isn't the closest node to xi
                    if( xi >= (xi_nodes[k] + xi_nodes[k+1])/2 ){
                        lskip *= (xi - xi_nodes[k]);
                    } else {
                        break; // k is the closest node to xi
                    }
                }
                for(int i = k + 1; i < Pn + 1; ++i){
                    lskip *= (xi - xi_nodes[i]);
                }
                T lprod = lskip * (xi - xi_nodes[k]);

                // calculate the sum of inverse differences
                // neglecting the skipped node
                // And calculate Nj in the same loops
                T s = 0.0;
                int j;
                for(j = 0; j < k; ++j){
                    T inv_diff = 1.0 / (xi - xi_nodes[j]);
                    s += inv_diff;
                    Nj[j] = lprod * inv_diff * wj[j];
                }
                Nj[k] = lskip * wj[k];
                for(++j; j < Pn + 1; ++j){
                    T inv_diff = 1.0 / (xi - xi_nodes[j]);
                    s += inv_diff;
                    Nj[j] = lprod * inv_diff * wj[j];
                }

                // run-time precompute the derivative of the l-product 
                T lprime = lprod * s + lskip;

                // evaluate the derivatives
                for(j = 0; j < k; ++j){
                    // quotient rule
                    dNj[j] = (lprime * wj[j] - Nj[j]) / (xi - xi_nodes[j]);
                }
                dNj[k] = s * Nj[k];
                for(++j; j < Pn + 1; ++j){
                    // quotient rule
                    dNj[j] = (lprime * wj[j] - Nj[j]) / (xi - xi_nodes[j]);
                }
            }
        }

        /**
         * @brief get the value, derivative, and second derivative 
         * of every interpolating polynomial at the given point
         * @param xi the point to evaluate at 
         * @param [out] Nj the basis function evaluations
         * @param [out] dNj the derivative evaluations
         * @param [out] d2Nj the second derivative evaluations
         */
        void d2_all(
            T xi,
            Tensor<T, Pn+1> &Nj,
            Tensor<T, Pn+1> &dNj,
            Tensor<T, Pn+1> &d2Nj
        ) const {
            deriv_all(xi, Nj, dNj);
            for(int j = 0; j < nbasis; ++j){
                d2Nj[j] = POLYNOMIAL::dNlagrange1d<T, Pn>(j, 2, xi);
            }
        }
    };
}
