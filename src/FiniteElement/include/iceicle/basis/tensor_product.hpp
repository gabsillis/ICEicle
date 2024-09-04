/**
 * @brief generalized interfaces for tensor products of basis functions
 * @author Gianni Absillis (gabsill@ncsu.edu)
 */
#pragma once

#include "Numtool/integer_utils.hpp"
#include "Numtool/fixed_size_tensor.hpp"
#include "Numtool/tmp_flow_control.hpp"
#include "Numtool/constexpr_math.hpp"
#include <Numtool/point.hpp>
#include <algorithm>
#include <cstddef>
#include <mdspan/mdspan.hpp>
#include <stdexcept>
#include <string>
#include <sstream>

namespace iceicle{

    /// @brief Cartesian Product of indicies 
    /// defined by the size of the indices in each dimension
    ///
    /// e.g CartesianIndexProduct{std::array{3, 2, 2}} represents 
    ///     {0, 0, 0}
    ///     {0, 0, 1}
    ///     {0, 1, 0}
    ///     ...
    ///     {2, 1, 1}
    ///
    /// @tparam index_type the data type of the indices 
    /// @tparam ndim the number of dimensions 
    /// @param sizes the extents of the indices in each dimension
    template<class index_type, std::size_t ndim>
    constexpr
    auto cartesian_index_product(std::array<index_type, ndim> sizes) noexcept 
    -> std::vector<std::array<index_type, ndim>>
    {
        std::vector<std::array<index_type, ndim>> ijk_poin{};
        std::array<index_type, ndim> ijk;
        std::ranges::fill(ijk, 0);
        ijk_poin.push_back(ijk);

        auto increment = [sizes](std::array<index_type, ndim> &ijk) {
            ijk[0]++;
            for(std::size_t idim = 0; idim < ndim - 1; ++idim) {
                if(ijk[idim] == sizes[idim]){
                    ijk[idim] = 0;
                    ijk[idim+1]++;
                } else {
                    break;
                }
            }

            if(ijk[ndim - 1] == sizes[ndim - 1]) 
                return false;
            else 
                return true;
        };

        while(increment(ijk)) ijk_poin.push_back(ijk);
        return ijk_poin;
    }

    /**
     * @brief represents a multi-index set for a Q-type tensor product
     * @tparam index_type the type used for indexing
     * @tparam ndim the number of dimensions 
     * @param size_1d the number of indices in a single dimension
     *
     * serves as a range of type std::array<index_type, ndim> 
     * so that all the multidimensional indices can be iterated over
     */
    template<class index_type, index_type ndim, index_type size_1d>
    struct QTypeIndexSet {

        // ============
        // = Typedefs =
        // ============

        using value_type = std::array<index_type, ndim>;
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using reference = value_type&;
        using const_reference = const value_type&;
        using pointer = value_type*;
        using const_pointer = const value_type*;
        using iterator = value_type*;
        using const_iterator = const value_type*;
        using reverse_iterator = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;

        // ================
        // = Data Members =
        // ================

        /// The total number of entries generated by the tensor product
        static constexpr size_type cardinality = MATH::power_T<size_1d, ndim>::value;

        /// Basis function indices by dimension for each node
        static constexpr std::array<std::array<index_type, ndim>, cardinality> ijk_poin = []() {
            std::array<std::array<index_type, ndim>, cardinality> ret{};

            NUMTOOL::TMP::constexpr_for_range<0, ndim>([&ret]<int idim>() {
                // number of times to repeat the loop over basis functions
                const int nrepeat = MATH::power_T<size_1d, idim>::value;
                // the size that one loop through the basis function indices gives
                const int cyclesize = MATH::power_T<size_1d, ndim - idim>::value;
                for (int irep = 0; irep < nrepeat; ++irep) {
                    NUMTOOL::TMP::constexpr_for_range<0, size_1d>(
                        [irep, &ret]<int ibasis>() {
                        const int nfill = NUMTOOL::TMP::pow(size_1d, ndim - idim - 1);

                        // offset for multiplying by this ibasis
                        const int start_offset = ibasis * nfill;

                        // multiply the next nfill by the current basis function
                        for (int ifill = 0; ifill < nfill; ++ifill) {
                            const int offset = irep * cyclesize + start_offset;
                            ret[offset + ifill][idim] = ibasis;
                        }
                    });
                }
            });

            return ret;
        }();

        // ==================
        // = Element Access =
        // ==================

        /// @brief get a reference to the multi-index at the given position
        /// @throws std::out_of_range if not in the cardinality of the multi-index set1
        inline constexpr 
        auto at(index_type pos) -> reference 
        {
            if(pos < 0 || pos >= cardinality) throw std::out_of_range("Index out of bounds.");
            return ijk_poin[pos];
        }

        /// @brief get a reference to the multi-index at the given position
        /// @throws std::out_of_range if not in the cardinality of the multi-index set1
        inline constexpr 
        auto at(index_type pos) const -> const_reference 
        {
            if(pos < 0 || pos >= cardinality) throw std::out_of_range("Index out of bounds.");
            return ijk_poin[pos];
        }

        /// @brief get a reference to the multi-index at the given position
        inline constexpr
        auto operator[] (index_type pos) noexcept -> reference 
        {
            return ijk_poin[pos];
        }

        /// @brief get a reference to the multi-index at the given position
        inline constexpr 
        auto operator[] (index_type pos) const noexcept -> const_reference 
        {
            return ijk_poin[pos];
        }

        /// @brief access the underlying array of multi-indices
        inline constexpr 
        auto data() noexcept -> pointer {
            return ijk_poin.data();
        }

        /// @brief access the underlying array of multi-indices
        inline constexpr 
        auto data() const noexcept -> const_pointer {
            return ijk_poin.data();
        }

        // =============
        // = Iterators =
        // =============
        
        /// @brief iterator to first multi-index
        inline constexpr 
        auto begin() noexcept -> iterator {
            return ijk_poin.begin();
        }

        /// @brief iterator to first multi-index
        inline constexpr 
        auto begin() const noexcept -> const_iterator {
            return ijk_poin.begin();
        }


        /// @brief iterator to first multi-index
        inline constexpr 
        auto cbegin() const noexcept -> const_iterator {
            return ijk_poin.begin();
        }

        /// @brief iterator to last multi-index
        inline constexpr 
        auto end() noexcept -> iterator {
            return ijk_poin.end();
        }

        /// @brief iterator to last multi-index
        inline constexpr 
        auto end() const noexcept -> const_iterator {
            return ijk_poin.end();
        }

        /// @brief iterator to last multi-index
        inline constexpr 
        auto cend() const noexcept -> const_iterator {
            return ijk_poin.end();
        }

        /// @brief iterator to first multi-index
        /// in the reverse sequence
        inline constexpr 
        auto rbegin() noexcept -> iterator {
            return ijk_poin.rbegin();
        }

        /// @brief iterator to first multi-index
        /// in the reverse sequence
        inline constexpr 
        auto rbegin() const noexcept -> const_iterator {
            return ijk_poin.rbegin();
        }


        /// @brief iterator to first multi-index
        /// in the reverse sequence
        inline constexpr 
        auto crbegin() const noexcept -> const_iterator {
            return ijk_poin.rbegin();
        }

        /// @brief iterator to last multi-index
        /// in the reverse sequence
        inline constexpr 
        auto rend() noexcept -> iterator {
            return ijk_poin.rend();
        }

        /// @brief iterator to last multi-index
        /// in the reverse sequence
        inline constexpr 
        auto rend() const noexcept -> const_iterator {
            return ijk_poin.rend();
        }

        /// @brief iterator to last multi-index
        /// in the reverse sequence
        inline constexpr 
        auto crend() const noexcept -> const_iterator {
            return ijk_poin.rend();
        }

        /// @brief checks if there are no elements in the multi-index set 
        [[nodiscard]] inline constexpr  
        auto empty() const noexcept -> bool { return begin() == end(); }

        /// @brief bet the number of elements in the multi-index set
        inline constexpr 
        auto size() const noexcept -> size_type { return cardinality; }
    };

    template<class Basis1D>
    concept basis_C0 = requires(const Basis1D &basis, Basis1D::value_type abscisse) {
        // The basis function must have a statically determinable number of basis functions
        {Basis1D::nbasis} -> std::convertible_to<int>;

        // The basis function set must be evaluatable
        {basis.eval_all(abscisse)} -> 
            std::convertible_to<
                NUMTOOL::TENSOR::FIXED_SIZE::Tensor<
                    typename Basis1D::value_type, 
                    Basis1D::nbasis
                >
            >;
    };

    template<class Basis1D>
    concept basis_C1 = basis_C0<Basis1D> && requires(
        const Basis1D &basis,
        Basis1D::value_type abscisse,
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<typename Basis1D::value_type, Basis1D::nbasis> &Nj,
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<typename Basis1D::value_type, Basis1D::nbasis> &dNj
    ) {
        // The basis function must have a statically determinable number of basis functions
        {Basis1D::nbasis} -> std::convertible_to<int>;

        // The basis function set must be evaluatable
        {basis.deriv_all(abscisse, Nj, dNj)};
    };

    template<class Basis1D>
    concept basis_C2 = basis_C0<Basis1D> && basis_C1<Basis1D> && requires(
        const Basis1D &basis,
        Basis1D::value_type abscisse,
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<typename Basis1D::value_type, Basis1D::nbasis> &Nj,
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<typename Basis1D::value_type, Basis1D::nbasis> &dNj,
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<typename Basis1D::value_type, Basis1D::nbasis> &d2Nj
    ) {
        // The basis function must have a statically determinable number of basis functions
        {Basis1D::nbasis} -> std::convertible_to<int>;

        // The basis function set must be evaluatable
        {basis.d2_all(abscisse, Nj, dNj, d2Nj)};
    };

    /**
     * @brief a Q-Type tensor product 
     * The classical outer product tensor product 
     * @tparam T the floating point type
     * @tparam ndim the number of dimensions
     * @tparam nbasis_1d the number of basis functions for the 1D basis function
     */
    template<class T, int ndim, int nbasis_1d>
    struct QTypeProduct {

        // === Aliases ===
        using Point = MATH::GEOMETRY::Point<T, ndim>;
        template<typename T1, std::size_t... sizes>
        using Tensor = NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T1, sizes...>;

        // =============
        // = Constants =
        // =============

        /// The total number of entries generated by the tensor product
        static constexpr int nvalues = MATH::power_T<nbasis_1d, ndim>::value;

        /// Basis function indices by dimension for each node
        static constexpr std::array<std::array<int, ndim>, nvalues> ijk_poin = []() {
            std::array<std::array<int, ndim>, nvalues> ret{};

            NUMTOOL::TMP::constexpr_for_range<0, ndim>([&ret]<int idim>() {
                // number of times to repeat the loop over basis functions
                const int nrepeat = MATH::power_T<nbasis_1d, idim>::value;
                // the size that one loop through the basis function indices gives
                const int cyclesize = MATH::power_T<nbasis_1d, ndim - idim>::value;
                for (int irep = 0; irep < nrepeat; ++irep) {
                    NUMTOOL::TMP::constexpr_for_range<0, nbasis_1d>(
                        [irep, &ret]<int ibasis>() {
                        const int nfill = NUMTOOL::TMP::pow(nbasis_1d, ndim - idim - 1);

                        // offset for multiplying by this ibasis
                        const int start_offset = ibasis * nfill;

                        // multiply the next nfill by the current basis function
                        for (int ifill = 0; ifill < nfill; ++ifill) {
                            const int offset = irep * cyclesize + start_offset;
                            ret[offset + ifill][idim] = ibasis;
                        }
                    });
                }
            });

            return ret;
        }();

        /// The distance in 1d indices between subsequent ijk (dimensional) indices 
        static constexpr NUMTOOL::TENSOR::FIXED_SIZE::Tensor<int, ndim> strides = []{
            NUMTOOL::TENSOR::FIXED_SIZE::Tensor<int, ndim> ret{};
            NUMTOOL::TMP::constexpr_for_range<0, ndim>([&]<int i>{
                ret[i] = MATH::power_T<nbasis_1d, ndim - i - 1>::value;
            });
            return ret;
        }();

        // ==================
        // = Multidim Index =
        // =     Utility    =
        // ==================

        /// @brief convert an ijk (dimensional) index to a 1d index
        static constexpr int convert_ijk(int ijk[ndim]) {
            int ret = 0;
            for(int idim = 0; idim < ndim; ++idim){
                ret += ijk[idim] * std::pow(nbasis_1d, ndim - idim - 1);
            }
            return ret;
        }


        /** @brief print the 1d lagrange basis function indices for each dimension for
        * each node */
        std::string print_ijk_poin() {
            using namespace std;
            std::ostringstream ijk_string;
            for (int inode = 0; inode < nvalues; ++inode) {
            ijk_string << "[";
            for (int idim = 0; idim < ndim; ++idim) {
                ijk_string << " " << ijk_poin[inode][idim];
            }
            ijk_string << " ]\n";
            }
            return ijk_string.str();
        }

        // ==================
        // = Tensor Product =
        // =   Evaluation   =
        // ==================


        /**
        * @brief fill the array with shape functions at the given point
        *
        * Requirements:
        *   - C0 basis function - the 0th derivatives are calculable in bulk
        *   - Value type of the basis function must match T
        *
        * @param [in] Basis1D the 1d basis functions
        * @param [in] xi the point in the reference domain to evaluate the basis at
        * @param [out] Bi the shape function evaluations
        */
        inline void fill_shp(
            const basis_C1 auto &basis_1d,
            const Point &xi,
            T *Bi
        ) const noexcept {
            using namespace NUMTOOL::TENSOR::FIXED_SIZE;

            if constexpr(ndim == 0){
                Bi[0] = 1;
            } else {
                // run-time precompute the lagrange polynomial evaluations 
                // for each coordinate
                Tensor<T, ndim, nbasis_1d> lagrange_evals{};
                for(int idim = 0; idim < ndim; ++idim){
                    lagrange_evals[idim] = basis_1d.eval_all(xi[idim]);
                }

                // for the first dimension (fencepost)
                NUMTOOL::TMP::constexpr_for_range<0, nbasis_1d>(
                    [&]<int ibasis>(T xi_dim) {
                        static constexpr int nfill = MATH::power_T<nbasis_1d, ndim - 1>::value;
                        T Bi_idim = lagrange_evals[0][ibasis];
                        std::fill_n(Bi + nfill * ibasis, nfill, Bi_idim);
                    },
                    xi[0]
                );

                for (int idim = 1; idim < ndim; ++idim) {
                    T xi_dim = xi[idim];

                    // number of times to repeat the loop over basis functions
                    int nrepeat = std::pow(nbasis_1d, idim);
                    // the size that one loop through the basis function indices gives
                    const int cyclesize = std::pow(nbasis_1d, ndim - idim);
                    for (int irep = 0; irep < nrepeat; ++irep) {
                        NUMTOOL::TMP::constexpr_for_range<0, nbasis_1d>(
                            [&]<int ibasis>(int idim, T xi_dim, T *Bi) {
                            // evaluate the 1d basis function at the idimth coordinate
                            T Bi_idim = lagrange_evals[idim][ibasis];
                            const int nfill = std::pow(nbasis_1d, ndim - idim - 1);

                            // offset for multiplying by this ibasis
                            const int start_offset = ibasis * nfill;

                            // multiply the next nfill by the current basis function
                            for (int ifill = 0; ifill < nfill; ++ifill) {
                                    Bi[start_offset + ifill] *= Bi_idim;
                                }
                            },
                            idim, xi_dim, Bi + irep * cyclesize
                        );
                    }
                }
            }
        }

        /**
        * @brief fill the given 2d array with the derivatives of each basis function 
        * @param [in] basis_1d the 1d basis functions
        * @param [in] xi the point in the reference domain to evaluate the derivative at 
        * @param [out] dBidxj the derivatives of the basis functions 
        */
        void fill_deriv(
            const basis_C1 auto &basis_1d,
            const Point &xi,
            Tensor<T, nvalues, ndim> &dBidxj
        ) const noexcept {
            using namespace NUMTOOL::TENSOR::FIXED_SIZE;

            if constexpr(ndim == 0){
                dBidxj.data()[0] = 0.0;
            } else {
                // run-time precompute the lagrange polynomial evaluations 
                // and derivatives for each coordinate
                Tensor<T, ndim, nbasis_1d> lagrange_evals{};
                Tensor<T, ndim, nbasis_1d> lagrange_derivs{};
                for(int idim = 0; idim < ndim; ++idim){
                basis_1d.deriv_all(
                    xi[idim], lagrange_evals[idim], lagrange_derivs[idim]);
                }
                // fencepost the loop at idim = 0
                NUMTOOL::TMP::constexpr_for_range<0, nbasis_1d>(
                [&]<int ibasis>(const Point &xi) {
                    static constexpr int nfill = MATH::power_T<nbasis_1d, ndim - 1>::value;
                    T Bi_idim = lagrange_evals[0][ibasis];
                    T dBi_idim = lagrange_derivs[0][ibasis];
                    for(int ifill = 0; ifill < nfill; ++ifill){
                        dBidxj[nfill * ibasis + ifill][0] = dBi_idim;
                        for(int jdim = 1; jdim < ndim; ++jdim){
                            dBidxj[nfill * ibasis + ifill][jdim] = Bi_idim;
                        }
                    }
                },
                xi);
                
                NUMTOOL::TMP::constexpr_for_range<1, ndim>(
                    [&]<int idim>(const Point &xi){
                        // number of times to repeat the loop over basis functions
                        const int nrepeat = std::pow(nbasis_1d, idim);
                        // the size that one loop through the basis function indices gives 
                        const int cyclesize = std::pow(nbasis_1d, ndim - idim);

                        for(int irep = 0; irep < nrepeat; ++irep) {
                            for(int ibasis = 0; ibasis < nbasis_1d; ++ibasis){
                                T dBi_idim = lagrange_derivs[idim][ibasis];
                                const int nfill = std::pow(nbasis_1d, ndim - idim - 1);

                                // offset for multiplying by this ibasis
                                const int start_offset = ibasis * nfill;

                                // multiply the next nfill by the current basis function
                                for (int ifill = 0; ifill < nfill; ++ifill) {
                                    const int offset = irep * cyclesize + start_offset;
                                    NUMTOOL::TMP::constexpr_for_range<0, ndim>([&]<int jdim>(){
                                        if constexpr(jdim == idim){
                                            dBidxj[offset + ifill][jdim] *= dBi_idim;
                                        } else {
                                            dBidxj[offset + ifill][jdim] *= lagrange_evals[idim][ibasis];
                                        }
                                    });
                                }
                            }
                        }
                    },
                    xi
                );
            }
        }

        /**
        * @brief fill the provided 1d array with the hessian of each basis function 
        * \frac{ d^2Bi }{ dx_j dx_k}
        * in C row major order this corresponds to [i][j][k]
        * and is symmetric in the last two indices 
        *
        * @param basis_1d [in] the basis function to get the tensor product for
        * @param [in] xi the point in the reference domain to evaluate at 
        * @param [out] nodal_hessian_data the pointer to fill with the hessian 
        * must be preallocated to n_nodes() * ndim * ndim 
        * (this get's zero'd out by this function before use)
        *
        * @return an mdspan view of nodal_hessian_data
        */
        auto fill_hess(
            const basis_C2 auto &basis_1d,
            const T *xi, T *nodal_hessian_data) const noexcept{
            using namespace std::experimental;
            using namespace NUMTOOL::TENSOR::FIXED_SIZE;

            // view the hessian output array by an mdspan 
            mdspan hess{nodal_hessian_data, extents{nvalues, ndim, ndim}};


            // special case fill
            if constexpr(ndim == 0){
                nodal_hessian_data[0] = 0.0;
                return hess;
            }
            // fill with ones for multiplicative identity
            std::fill_n(nodal_hessian_data, nvalues*ndim*ndim, 1.0);
            
            // run-time precompute the lagrange polynomial evaluations 
            // and derivatives for each coordinate
            Tensor<T, ndim, nbasis_1d> lagrange_evals{};
            Tensor<T, ndim, nbasis_1d> lagrange_derivs{};
            Tensor<T, ndim, nbasis_1d> lagrange_d2s{};
            for(int idim = 0; idim < ndim; ++idim){
                basis_1d.d2_all(
                    xi[idim], lagrange_evals[idim],
                    lagrange_derivs[idim], lagrange_d2s[idim]);
            }

            // TODO: optimize using the fill strategy like fill_deriv 
            for(int ibasis = 0; ibasis < nvalues; ++ibasis){

                for(int ideriv = 0; ideriv < ndim; ++ideriv){
                    for(int jderiv = ideriv; jderiv < ndim; ++jderiv){
                        // handle diagonal terms 
                        if(ideriv == jderiv){
                        for(int idim = 0; idim < ndim; ++idim){
                            if(idim == ideriv){
                            hess[ibasis, ideriv, jderiv] *= lagrange_d2s[idim][ijk_poin[ibasis][idim]];
                            } else {
                            hess[ibasis, ideriv, jderiv] *= lagrange_evals[idim][ijk_poin[ibasis][idim]];
                            }
                        }
                        } else {
                        // loop over the 1d basis functions for each dimension
                        for(int idim = 0; idim < ndim; ++idim){
                            if(idim == ideriv){
                            hess[ibasis, ideriv, jderiv] *= 
                                lagrange_derivs[ideriv][ijk_poin[ibasis][idim]];
                            } else if(idim == jderiv){
                            hess[ibasis, ideriv, jderiv] *= 
                                lagrange_derivs[jderiv][ijk_poin[ibasis][idim]];
                            } else {
                            // not a derivative so just the 1d function 
                            hess[ibasis, ideriv, jderiv] *=
                                lagrange_evals[idim][ijk_poin[ibasis][idim]];
                            }
                        }
                        }
                    }
                }

                // copy symmetric part 
                for(int ideriv = 0; ideriv < ndim; ++ideriv){
                    for(int jderiv = 0; jderiv < ideriv; ++jderiv){
                        hess[ibasis, ideriv, jderiv] = hess[ibasis, jderiv, ideriv];
                    }
                }
            }

            return hess;
        }
    };

}