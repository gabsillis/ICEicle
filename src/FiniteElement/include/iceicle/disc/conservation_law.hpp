#pragma once

#include "Numtool/fixed_size_tensor.hpp"
#include "Numtool/point.hpp"
#include "iceicle/element/finite_element.hpp"
#include "iceicle/fd_utils.hpp"
#include "iceicle/fe_function/fespan.hpp"
#include "iceicle/geometry/face.hpp"
#include "iceicle/mesh/mesh.hpp"
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>
namespace iceicle {

    /// @brief a Physical flux given a vector valued state u 
    /// and the gradient of u 
    /// returns a vector valued flux for each dimension F(u, gradu)
    template<class FluxT>
    concept physical_flux = 
    requires(
        FluxT flux,
        std::array<typename FluxT::value_type, FluxT::nv_comp> u,
        std::mdspan<typename FluxT::value_type, std::extents<std::size_t, std::dynamic_extent, std::dynamic_extent>> gradu
    ) {
        { flux(u, gradu) } -> std::same_as<NUMTOOL::TENSOR::FIXED_SIZE::Tensor<typename FluxT::value_type, FluxT::nv_comp, FluxT::ndim>> ;
    };

    /// @brief Implementation of a numerical flux for convective fluxes 
    ///
    /// given a state uL and uR on either side of an interface and the unit normal vector of the interface 
    /// return the flux
    template<class FluxT>
    concept convective_numerical_flux = requires(
        FluxT flux,
        std::array<typename FluxT::value_type, FluxT::nv_comp> uL,
        std::array<typename FluxT::value_type, FluxT::nv_comp> uR,
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor< typename FluxT::value_type, FluxT::ndim > unit_normal
    ) {
        { flux(uL, uR, unit_normal) } -> std::same_as<std::array<typename FluxT::value_type, FluxT::nv_comp>>;
    };

    /// @brief the diffusive flux normal to the interface 
    ///
    /// given a single valued solution at an interface and the single valued gradient 
    /// compute the flux function for diffusion operators in the normal direction
    ///
    /// NOTE: this will be evaluated separately from ConvectiveNumericalFlux so do not include that here
    template<class FluxT>
    concept diffusion_flux = 
    requires(
        FluxT flux,
        std::array<typename FluxT::value_type, FluxT::nv_comp> u,
        std::mdspan<typename FluxT::value_type, std::extents<std::size_t, std::dynamic_extent, std::dynamic_extent>> gradu,
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor< typename FluxT::value_type, FluxT::ndim > unit_normal
    ) {
        { flux(u, gradu, unit_normal) } -> std::same_as<std::array<typename FluxT::value_type, FluxT::nv_comp>>;
    };

    /// @brief additional boundary conditions that are implemented at the PDE level
    /// This includes things like characteristic and wall boundary conditions 
    template<class FluxT>
    concept implements_bcs =
    requires(
        FluxT flux,
        std::array<typename FluxT::value_type, FluxT::nv_comp> &uL,
        std::mdspan<typename FluxT::value_type, std::extents<std::size_t, std::dynamic_extent, std::dynamic_extent>> graduL,
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor< typename FluxT::value_type, FluxT::ndim > unit_normal,
        BOUNDARY_CONDITIONS bctype,
        int bcflag
    ) {
        { flux.apply_bc(uL, graduL, unit_normal, bctype, bcflag)} -> std::same_as<std::pair< 
                std::array<typename FluxT::value_type, FluxT::nv_comp>, 
                NUMTOOL::TENSOR::FIXED_SIZE::Tensor< typename FluxT::value_type, 
                    FluxT::neq, FluxT::ndim>
            >>;
                
    };

    /// @brief Diffusion fluxes can explicitly define the homogeineity tensor given a state u 
    /// This can be used for interface correction
    template< class FluxT>
    concept computes_homogeneity_tensor = requires(
        FluxT flux,
        std::array<typename FluxT::value_type, FluxT::nv_comp> u
    ) {
        { flux.homogeneity_tensor(u) } -> std::same_as< 
                NUMTOOL::TENSOR::FIXED_SIZE::Tensor<typename FluxT::value_type, FluxT::nv_comp, FluxT::ndim, FluxT::nv_comp, FluxT::ndim>>;
    };

    template<
        typename T,
        int ndim,
        physical_flux PFlux,
        convective_numerical_flux CFlux,
        diffusion_flux DiffusiveFlux,
        class ST_Info = std::false_type
    >
    class ConservationLawDDG {

        public:
        PFlux phys_flux;
        CFlux conv_nflux;
        DiffusiveFlux diff_flux;

        // compile time correctness checks
        static_assert(std::same_as<T, typename PFlux::value_type>, "Value type mismatch with physical flux");
        static_assert(std::same_as<T, typename CFlux::value_type>, "Value type mismatch with convective flux");
        static_assert(std::same_as<T, typename DiffusiveFlux::value_type>, "Value type mismatch with diffusive flux");
        static_assert(ndim == PFlux::ndim, "ndim mismatch with physical flux");
        static_assert(ndim == CFlux::ndim, "ndim mismatch with convective flux");
        static_assert(ndim == DiffusiveFlux::ndim, "ndim mismatch with diffusive flux");

        // ============
        // = Typedefs =
        // ============

        using value_type = T;

        // =============
        // = Constants =
        // =============

        /// @brief access the number of dimensions through a public interface
        static constexpr int dimensionality = ndim;

        /// @brief the number of vector components
        static constexpr std::size_t nv_comp = PFlux::nv_comp;
        static constexpr std::size_t dnv_comp = PFlux::nv_comp;

        // ==================
        // = Public Members =
        // ==================

        /// @brief switch to use the interior penalty method instead of ddg 
        bool interior_penalty = false;

        /// @brief IC multiplier to get DDGIC
        /// see Danis Yan 2023 Journal of Scientific Computing
        /// DDGIC (sigma = 1)
        /// Default: Standard DDG (sigma = 0)
        T sigma_ic = 0.0;

        /// @brief dirichlet value for each bcflag index
        /// as a function callback 
        /// This function will take the physical domain point (size = ndim)
        /// and output neq values in the second argument
        std::vector< std::function<void(const T*, T*)> > dirichlet_callbacks;

        /// @brief neumann value for each bcflag index
        /// as a function callback 
        /// This function will take the physical domain point (size = ndim)
        /// and output neq values in the second argument
        std::vector< std::function<void(const T*, T*)> > neumann_callbacks;

        /// @brief user defined source term as a function callback 
        /// Takes the position in the domain (size = ndim) 
        /// and outputs the source (size = neq) in the second argument
        std::optional< std::function<void(const T*, T*)> > user_source = std::nullopt;

        /// @brief utility for SPACETIME_PAST boundary condition
        ST_Info spacetime_info;

        /// @brief human readable names for each vector component of the variables
        std::vector<std::string> field_names;

        /// @brief human readable names for each vector component of the residuals
        std::vector<std::string> residual_names;

        // ===============
        // = Constructor =
        // ===============
        //
        /** @brief construct from 
         * and take ownership of the fluxes and spacetime_info
         *
         * @param physical_flux the actual discretization flux 
         * F in grad(F) = 0 
         * or for method of lines:
         * du/dt + grad(F) = 0
         *
         * @param convective_numflux the numerical flux for the 
         * convective portion (typically a Riemann Solver or Upwinding method)
         *
         * @param diffusive_numflux the numerical flux for the diffusive portion
         *
         * @param spacetime_info is the class that defines the SPACETIME_PAST connection
         * see iceicle::SpacetimeConnection
         */
        constexpr ConservationLawDDG(
            PFlux&& physical_flux,
            CFlux&& convective_numflux,
            DiffusiveFlux&& diffusive_flux,
            ST_Info&& spacetime_info
        ) noexcept : phys_flux{physical_flux}, conv_nflux{convective_numflux}, 
            diff_flux{diffusive_flux}, spacetime_info{spacetime_info} {}


        /** @brief construct from 
         * and take ownership of the fluxes
         *
         * @param physical_flux the actual discretization flux 
         * F in grad(F) = 0 
         * or for method of lines:
         * du/dt + grad(F) = 0
         *
         * @param convective_numflux the numerical flux for the 
         * convective portion (typically a Riemann Solver or Upwinding method)
         *
         * @param diffusive_numflux the numerical flux for the diffusive portion
         */
        constexpr ConservationLawDDG(
            PFlux&& physical_flux,
            CFlux&& convective_numflux,
            DiffusiveFlux&& diffusive_flux
        ) noexcept : phys_flux{physical_flux}, conv_nflux{convective_numflux}, 
            diff_flux{diffusive_flux} {}


        // ============================
        // = Discretization Interface =
        // ============================

        /**
         * @brief get the timestep from cfl 
         * this takes it from the physical flux
         * often this will require data to be set from the domain and boundary integrals 
         * such as wavespeeds, which will arise naturally during residual computation
         * (WARNING: except for the very first iteration)
         *
         * @param cfl the cfl condition 
         * @param reference_length the size to use for the length of the cfl condition 
         * @return the timestep based on the cfl condition
         */
        T dt_from_cfl(T cfl, T reference_length){
            return phys_flux.dt_from_cfl(cfl, reference_length);
        }

        // =============
        // = Integrals =
        // =============

        template<class IDX>
        auto domain_integral(
            const FiniteElement<T, IDX, ndim> &el,
            elspan auto unkel,
            elspan auto res
        ) const -> void {
            static constexpr int neq = decltype(unkel)::static_extent();
            static_assert(neq == PFlux::nv_comp, "Number of equations must match.");
            using namespace NUMTOOL::TENSOR::FIXED_SIZE;

            // basis function scratch space
            std::vector<T> dbdx_data(el.nbasis() * ndim);

            // solution scratch space
            std::array<T, neq> u;
            std::array<T, neq * ndim> gradu_data;

            // loop over the quadrature points
            for(int iqp = 0; iqp < el.nQP(); ++iqp){
                const QuadraturePoint<T, ndim> &quadpt = el.getQP(iqp);

                // calculate the jacobian determinant 
                auto J = el.jacobian(quadpt.abscisse);
                T detJ = NUMTOOL::TENSOR::FIXED_SIZE::determinant(J);
                // prevent duplicate contribution of overlapping range in transformation
                // this occurs in concave elements
                detJ = std::max((T) 0.0, detJ); 

                // get the basis functions and gradients in the physical domain
                auto bi = el.eval_basis_qp(iqp);
                auto gradxBi = el.eval_phys_grad_basis(quadpt.abscisse, J,
                        el.eval_grad_basis_qp(iqp), dbdx_data.data());

                // construct the solution U at the quadrature point 
                std::ranges::fill(u, 0.0);
                for(int ibasis = 0; ibasis < el.nbasis(); ++ibasis){
                    for(int ieq = 0; ieq < neq; ++ieq){
                        u[ieq] += unkel[ibasis, ieq] * bi[ibasis];
                    }
                }

                // construct the gradient of u 
                auto gradu = unkel.contract_mdspan(gradxBi, gradu_data.data());

                // compute the flux  and scatter to the residual
                Tensor<T, neq, ndim> flux = phys_flux(u, gradu);

                // loop over the test functions and construct the residual
                for(int itest = 0; itest < el.nbasis(); ++itest){
                    for(int ieq = 0; ieq < neq; ++ieq){
                        for(int jdim = 0; jdim < ndim; ++jdim){
                            res[itest, ieq] += flux[ieq][jdim] * gradxBi[itest, jdim] * detJ * quadpt.weight;
                        }
                    }
                }

                // if a source term has been defined, add it in
                if(user_source){
                    auto phys_pt = el.transform(quadpt.abscisse);

                    auto& source_fcn = user_source.value();
                    std::array<T, neq> source;
                    source_fcn(phys_pt.data(), source.data());
                    for(int itest = 0; itest < el.nbasis(); ++itest){
                        for(int ieq = 0; ieq < neq; ++ieq) {
                            res[itest, ieq] -= source[ieq] * bi[itest] * detJ * quadpt.weight;
                        }
                    }
                }
            }
        }

        template<class IDX>
        auto domain_integral_jacobian(
            const FiniteElement<T, IDX, ndim>& el,
            elspan auto unkel,
            linalg::out_matrix auto dfdu
        ) {
            static constexpr int neq = decltype(unkel)::static_extent();
            static_assert(neq == PFlux::nv_comp, "Number of equations must match.");
            using namespace NUMTOOL::TENSOR::FIXED_SIZE;
            // basis function scratch space
            std::vector<T> dbdx_data(el.nbasis() * ndim);

            // solution scratch space
            std::array<T, neq> u;
            std::array<T, neq * ndim> gradu_data;

            // get the degree of freedom layout for the element 
            auto el_layout = unkel.get_layout();

            // loop over the quadrature points
            for(int iqp = 0; iqp < el.nQP(); ++iqp){
                const QuadraturePoint<T, ndim> &quadpt = el.getQP(iqp);

                // calculate the jacobian determinant 
                auto J = el.jacobian(quadpt.abscisse);
                T detJ = NUMTOOL::TENSOR::FIXED_SIZE::determinant(J);
                // prevent duplicate contribution of overlapping range in transformation
                // this occurs in concave elements
                detJ = std::max((T) 0.0, detJ); 

                // get the basis functions and gradients in the physical domain
                auto bi = el.eval_basis_qp(iqp);
                auto gradxBi = el.eval_phys_grad_basis(quadpt.abscisse, J,
                        el.eval_grad_basis_qp(iqp), dbdx_data.data());

                // construct the solution U at the quadrature point 
                std::ranges::fill(u, 0.0);
                for(int ibasis = 0; ibasis < el.nbasis(); ++ibasis){
                    for(int ieq = 0; ieq < neq; ++ieq){
                        u[ieq] += unkel[ibasis, ieq] * bi[ibasis];
                    }
                }

                // construct the gradient of u 
                auto gradu = unkel.contract_mdspan(gradxBi, gradu_data.data());

                // compute the jacobian of the physical flux 
                // with respect to field variables
                // and gradients
                Tensor<T, neq, ndim> flux = phys_flux(u, gradu);
                T epsilon = scale_fd_epsilon(
                    std::sqrt(std::numeric_limits<T>::epsilon()),
                    frobenius(flux)
                );

                Tensor<T, neq, ndim, neq> dflux_du;
                Tensor<T, neq, ndim, neq, ndim> dflux_dgradu;
                for(int jeq = 0; jeq < neq; ++jeq){

                    // peturb u
                    T u_old = u[jeq];
                    u[jeq] += epsilon;
                    Tensor<T, neq, ndim> flux_p = phys_flux(u, gradu);
                    Tensor<T, neq, ndim> diff = (flux_p - flux) / epsilon;
                    for(int ieq = 0; ieq < neq; ++ieq){
                        for(int idim = 0; idim < ndim; ++idim){
                            dflux_du[ieq, idim, jeq] =
                                diff[ieq, idim];
                        }
                    }
                    u[jeq] = u_old;

                    // peturb grad u
                    for(int jdim = 0; jdim < ndim; ++jdim){
                        T gradu_old = gradu[jeq, jdim];
                        gradu[jeq, jdim] += epsilon;
                    
                        Tensor<T, neq, ndim> flux_p = phys_flux(u, gradu);
                        Tensor<T, neq, ndim> diff = (flux_p - flux) / epsilon;
                        for(int ieq = 0; ieq < neq; ++ieq){
                            for(int idim = 0; idim < ndim; ++idim){
                                dflux_dgradu[ieq, idim, jeq, jdim] =
                                    diff[ieq, idim];
                            }
                        }
                        gradu[jeq, jdim] = gradu_old;
                    }
                }
                
                // loop over the test functions and construct the jacobian
                for(int itest = 0; itest < el.nbasis(); ++itest){
                    for(int ieq = 0; ieq < neq; ++ieq){
                        for(int idim = 0; idim < ndim; ++idim){
                            for(int jdof = 0; jdof < el.nbasis(); ++jdof){
                                for(int jeq = 0; jeq < neq; ++jeq){
                                    // one-dimensional jacobian indices
                                    auto ijac = el_layout[itest, ieq];
                                    auto jjac = el_layout[jdof, jeq];
                                    dfdu[ijac, jjac] += 
                                        dflux_du[ieq, idim, jeq] * bi[jdof] 
                                        * gradxBi[itest, idim] * detJ * quadpt.weight;
                                    for(int jdim = 0; jdim < ndim; ++jdim){
                                        dfdu[ijac, jjac] += 
                                            dflux_dgradu[ieq, idim, jeq, jdim] * gradxBi[jdof, jdim] 
                                            * gradxBi[itest, idim] * detJ * quadpt.weight;
                                    }
                                }
                            }
                        }
                    }
                }

                // source has no dependence on U 
                // therefore no contribution to 
                // jacobian
            }
        }

        template<class IDX, class ULayoutPolicy, class UAccessorPolicy, class ResLayoutPolicy>
        void trace_integral(
            const TraceSpace<T, IDX, ndim> &trace,
            NodeArray<T, ndim> &coord,
            dofspan<T, ULayoutPolicy, UAccessorPolicy> unkelL,
            dofspan<T, ULayoutPolicy, UAccessorPolicy> unkelR,
            dofspan<T, ResLayoutPolicy> resL,
            dofspan<T, ResLayoutPolicy> resR
        ) const requires ( 
            elspan<decltype(unkelL)> && 
            elspan<decltype(unkelR)> && 
            elspan<decltype(resL)> && 
            elspan<decltype(resL)>
        ) {
            static constexpr int neq = nv_comp;
            static_assert(neq == decltype(unkelL)::static_extent(), "Number of equations must match.");
            static_assert(neq == decltype(unkelR)::static_extent(), "Number of equations must match.");
            using namespace NUMTOOL::TENSOR::FIXED_SIZE;
            using FiniteElement = FiniteElement<T, IDX, ndim>;

            // calculate the centroids of the left and right elements
            // in the physical domain
            const FiniteElement &elL = trace.elL;
            const FiniteElement &elR = trace.elR;
            auto centroidL = elL.centroid();
            auto centroidR = elR.centroid();

            // Basis function scratch space 
            PhysDomainEvalStorage storageL{elL};
            PhysDomainEvalStorage storageR{elR};

            // solution scratch space 
            std::array<T, neq> uL;
            std::array<T, neq> uR;
            std::array<T, neq * ndim> graduL_data;
            std::array<T, neq * ndim> graduR_data;
            std::array<T, neq * ndim> grad_ddg_data;
            std::array<T, neq * ndim * ndim> hessuL_data;
            std::array<T, neq * ndim * ndim> hessuR_data;


            // loop over the quadrature points 
            for(int iqp = 0; iqp < trace.nQP(); ++iqp){
                const QuadraturePoint<T, ndim - 1> &quadpt = trace.getQP(iqp);

                // calculate the riemannian metric tensor root
                auto Jfac = trace.face->Jacobian(coord, quadpt.abscisse);
                T sqrtg = trace.face->rootRiemannMetric(Jfac, quadpt.abscisse);

                // calculate the normal vector 
                auto normal = calc_ortho(Jfac);
                auto unit_normal = normalize(normal);

                // get the basis functions, derivatives, and hessians
                // (derivatives are wrt the physical domain)
                auto biL = trace.qp_evals_l[iqp].bi_span;
                auto biR = trace.qp_evals_r[iqp].bi_span;
                auto xiL = trace.transform_xiL(quadpt.abscisse);
                auto xiR = trace.transform_xiR(quadpt.abscisse);
                PhysDomainEval evalL{storageL, elL, xiL, trace.qp_evals_l[iqp]};
                PhysDomainEval evalR{storageR, elR, xiR, trace.qp_evals_r[iqp]};

                // construct the solution on the left and right
                std::ranges::fill(uL, 0.0);
                std::ranges::fill(uR, 0.0);
                for(int ieq = 0; ieq < neq; ++ieq){
                    for(int ibasis = 0; ibasis < elL.nbasis(); ++ibasis)
                        { uL[ieq] += unkelL[ibasis, ieq] * biL[ibasis]; }
                    for(int ibasis = 0; ibasis < elR.nbasis(); ++ibasis)
                        { uR[ieq] += unkelR[ibasis, ieq] * biR[ibasis]; }
                }

                // get the solution gradient and hessians
                auto graduL = unkelL.contract_mdspan(evalL.phys_grad_basis, graduL_data.data());
                auto graduR = unkelR.contract_mdspan(evalR.phys_grad_basis, graduR_data.data());
                auto hessuL = unkelL.contract_mdspan(evalL.phys_hess_basis, hessuL_data.data());
                auto hessuR = unkelR.contract_mdspan(evalR.phys_hess_basis, hessuR_data.data());

                // compute convective fluxes
                std::array<T, neq> fadvn = conv_nflux(uL, uR, unit_normal);

                // compute a single valued gradient using DDG or IP 

                // calculate the DDG distance
                MATH::GEOMETRY::Point<T, ndim> phys_pt;
                trace.face->transform(quadpt.abscisse, coord, phys_pt);
                T h_ddg = 0;
                for(int idim = 0; idim < ndim; ++idim){
                    h_ddg += unit_normal[idim] * (
                        (phys_pt[idim] - centroidL[idim])
                        + (centroidR[idim] - phys_pt[idim])
                    );
                }
                h_ddg = std::copysign(std::max(std::abs(h_ddg), std::numeric_limits<T>::epsilon()), h_ddg);
                
                int order = std::max(
                    elL.basis->getPolynomialOrder(),
                    elR.basis->getPolynomialOrder()
                );
                // Danis and Yan reccomended for NS
                T beta0 = std::pow(order + 1, 2);
                T beta1 = 1 / std::max((T) (2 * order * (order + 1)), 1.0);

                // switch to interior penalty if set
                if(interior_penalty) beta1 = 0.0;

                std::mdspan<T, std::extents<int, neq, ndim>> grad_ddg{grad_ddg_data.data()};
                for(int ieq = 0; ieq < neq; ++ieq){
                    // construct the DDG derivatives
                    T jumpu = uR[ieq] - uL[ieq];
                    for(int idim = 0; idim < ndim; ++idim){
                        grad_ddg[ieq, idim] = beta0 * jumpu / h_ddg * unit_normal[idim]
                            + 0.5 * (graduL[ieq, idim] + graduR[ieq, idim]);
                        T hessTerm = 0;
                        for(int jdim = 0; jdim < ndim; ++jdim){
                            hessTerm += (hessuR[ieq, jdim, idim] - hessuL[ieq, jdim, idim])
                                * unit_normal[jdim];
                        }
                        grad_ddg[ieq, idim] += beta1 * h_ddg * hessTerm;
                    }
                }

                // construct the viscous fluxes 
                std::array<T, neq> uavg;
                for(int ieq = 0; ieq < neq; ++ieq) uavg[ieq] = 0.5 * (uL[ieq] + uR[ieq]);

                std::array<T, neq> fviscn = diff_flux(uavg, grad_ddg, unit_normal);

                // scale by weight and face metric tensor
                for(int ieq = 0; ieq < neq; ++ieq){
                    fadvn[ieq] *= quadpt.weight * sqrtg;
                    fviscn[ieq] *= quadpt.weight * sqrtg;
                }

                // scatter contribution 
                for(int itest = 0; itest < elL.nbasis(); ++itest){
                    for(int ieq = 0; ieq < neq; ++ieq){
                        resL[itest, ieq] += (fviscn[ieq] - fadvn[ieq]) * biL[itest];
                    }
                }
                for(int itest = 0; itest < elR.nbasis(); ++itest){
                    for(int ieq = 0; ieq < neq; ++ieq){
                        resR[itest, ieq] -= (fviscn[ieq] - fadvn[ieq]) * biR[itest];
                    }
                }

                // if applicable: apply the interface correction 
                if constexpr (computes_homogeneity_tensor<DiffusiveFlux>) {
                    auto gradBiL = evalL.phys_grad_basis;
                    auto gradBiR = evalR.phys_grad_basis;
                    if(sigma_ic != 0.0){
                        auto Gtensor = diff_flux.homogeneity_tensor(uavg);

                        for(int itest = 0; itest < elL.nbasis(); ++itest){

                            for(int ieq = 0; ieq < neq; ++ieq){
                                for(int kdim = 0; kdim < ndim; ++kdim){
                                    for(int req = 0; req < neq; ++req){
                                        T jumpu_r = uR[req] - uL[req];
                                        for(int sdim = 0; sdim < ndim; ++sdim){
                                            
                                            T ic_contrib = 
                                                sigma_ic * Gtensor[ieq][kdim][req][sdim] * unit_normal[kdim]
                                                * jumpu_r * quadpt.weight * sqrtg;
                                            resL[itest, ieq] -= 
                                                ic_contrib
                                                * 0.5 * gradBiL[itest, sdim]; // 0.5 comes from average operator
                                            resR[itest, ieq] -= 
                                                ic_contrib
                                                * 0.5 * gradBiR[itest, sdim];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        /**
         * @brief calculate the weak form for a boundary condition 
         *        NOTE: Left is the interior element
         *
         * @tparam ULayoutPolicy the layout policy for the element data 
         * @tparam UAccessorPolicy the accessor policy for element data 
         * @tparam ResLayoutPolicy the layout policy for the residual data 
         *         (the accessor is the default for the residual)
         *
         * @param [in] trace the trace to integrate over 
         * @param [in] coord the global node coordinates array
         * @param [in] uL the interior element basis coefficients 
         * @param [in] uR is the same as uL unless this is a periodic boundary
         *                then this is the coefficients for the periodic element 
         * @param [out] resL the residual for the interior element 
         */
        template<class IDX, class ULayoutPolicy, class UAccessorPolicy, class ResLayoutPolicy>
        void boundaryIntegral(
            const TraceSpace<T, IDX, ndim> &trace,
            NodeArray<T, ndim> &coord,
            dofspan<T, ULayoutPolicy, UAccessorPolicy> unkelL,
            dofspan<T, ULayoutPolicy, UAccessorPolicy> unkelR,
            dofspan<T, ResLayoutPolicy> resL
        ) const requires(
            elspan<decltype(unkelL)> &&
            elspan<decltype(unkelR)> &&
            elspan<decltype(resL)> 
        ) {
            using namespace NUMTOOL::TENSOR::FIXED_SIZE;
            using FiniteElement = FiniteElement<T, IDX, ndim>;
            const FiniteElement &elL = trace.elL;

            static constexpr int neq = nv_comp;

            // Basis function scratch space 
            std::vector<T> gradbL_data(elL.nbasis() * ndim);

            // solution scratch space 
            std::array<T, neq> uL;
            std::array<T, neq * ndim> graduL_data;
            std::array<T, neq * ndim> graduR_data;
            std::array<T, neq * ndim> grad_ddg_data;
            std::array<T, neq * ndim * ndim> hessuL_data;
            std::array<T, neq * ndim * ndim> hessuR_data;
            auto centroidL = elL.centroid();

            // switch over special cases of boundary condition implementations 
            // that can have increased efficiency when coded separately
            switch(trace.face->bctype){
                case BOUNDARY_CONDITIONS::DIRICHLET: 
                {
                    // see Huang, Chen, Li, Yan 2016

                    // loop over quadrature points
                    for(int iqp = 0; iqp < trace.nQP(); ++iqp){
                        const QuadraturePoint<T, ndim - 1> &quadpt = trace.getQP(iqp);

                        // calculate the jacobian and riemannian metric root det
                        auto Jfac = trace.face->Jacobian(coord, quadpt.abscisse);
                        T sqrtg = trace.face->rootRiemannMetric(Jfac, quadpt.abscisse);

                        // calculate the normal vector 
                        auto normal = calc_ortho(Jfac);
                        auto unit_normal = normalize(normal);

                        // calculate the physical domain position
                        MATH::GEOMETRY::Point<T, ndim> phys_pt;
                        trace.face->transform(quadpt.abscisse, coord, phys_pt);

                        // get the function values
                        auto biL = trace.qp_evals_l[iqp].bi_span;

                        // get the gradients the physical domain
                        auto gradBiL = trace.eval_phys_grad_basis_l_qp(iqp, gradbL_data.data());
                        auto graduL = unkelL.contract_mdspan(gradBiL, graduL_data.data());

                        // construct the solution on the left and right
                        std::ranges::fill(uL, 0.0);
                        for(int ieq = 0; ieq < neq; ++ieq){
                            for(int ibasis = 0; ibasis < elL.nbasis(); ++ibasis)
                                { uL[ieq] += unkelL[ibasis, ieq] * biL[ibasis]; }
                        }

                        // Get the values at the boundary 
                        std::array<T, nv_comp> dirichlet_vals{};
                        dirichlet_callbacks[trace.face->bcflag](phys_pt.data(), dirichlet_vals.data());

                        // compute convective fluxes
                        std::array<T, neq> fadvn = conv_nflux(uL, dirichlet_vals, unit_normal);
                        // std::array<T, neq> fadvn = conv_nflux(dirichlet_vals, dirichlet_vals, unit_normal);

                        // calculate the DDG distance
                        T h_ddg = 0; // uses distance to quadpt on boundary face
                        for(int idim = 0; idim < ndim; ++idim){
                            h_ddg += std::abs(unit_normal[idim] * 
                                (phys_pt[idim] - centroidL[idim])
                            );
                        }
                        h_ddg = std::copysign(std::max(std::abs(h_ddg), std::numeric_limits<T>::epsilon()), h_ddg);

                        // construct the DDG derivatives
                        int order = 
                            elL.basis->getPolynomialOrder();
                        // Danis and Yan reccomended for NS
                        T beta0 = std::pow(order + 1, 2);
                        T beta1 = 1 / std::max((T) (2 * order * (order + 1)), 1.0);

                        std::mdspan<T, std::extents<int, neq, ndim>> grad_ddg{grad_ddg_data.data()};
                        for(int ieq = 0; ieq < neq; ++ieq){
                            // construct the DDG derivatives
                            T jumpu = dirichlet_vals[ieq] - uL[ieq];
                            for(int idim = 0; idim < ndim; ++idim){
                                grad_ddg[ieq, idim] = beta0 * jumpu / h_ddg * unit_normal[idim]
                                    + (graduL[ieq, idim]);
                            }
                        }

                        // construct the viscous fluxes 
                        std::array<T, neq> uavg;
                        for(int ieq = 0; ieq < neq; ++ieq) uavg[ieq] = 0.5 * (uL[ieq] + dirichlet_vals[ieq]);

                        std::array<T, neq> fviscn = diff_flux(uavg, grad_ddg, unit_normal);

                        // scale by weight and face metric tensor
                        for(int ieq = 0; ieq < neq; ++ieq){
                            fadvn[ieq] *= quadpt.weight * sqrtg;
                            fviscn[ieq] *= quadpt.weight * sqrtg;
                        }

                        // scatter contribution 
                        for(int itest = 0; itest < elL.nbasis(); ++itest){
                            for(int ieq = 0; ieq < neq; ++ieq){
                                resL[itest, ieq] += (fviscn[ieq] - fadvn[ieq]) * biL[itest];
                            }
                        }

                        // if applicable: apply the interface correction 
                        if constexpr (computes_homogeneity_tensor<DiffusiveFlux>) {
                            if(sigma_ic != 0.0){

                                auto Gtensor = diff_flux.homogeneity_tensor(uavg);

                                T average_gradv[ndim];
                                for(int itest = 0; itest < elL.nbasis(); ++itest){
                                    // get the average test function gradient
                                    for(int idim = 0; idim < ndim; ++idim){
                                        average_gradv[idim] = gradBiL[itest, idim];
                                    }

                                    for(int ieq = 0; ieq < neq; ++ieq){
                                        for(int kdim = 0; kdim < ndim; ++kdim){
                                            for(int req = 0; req < neq; ++req){
                                                T jumpu_r = dirichlet_vals[req] - uL[req];
                                                for(int sdim = 0; sdim < ndim; ++sdim){
                                                    
                                                    resL[itest, ieq] -= 
                                                        sigma_ic * Gtensor[ieq][kdim][req][sdim] * unit_normal[kdim] 
                                                        * average_gradv[sdim] * jumpu_r
                                                        * quadpt.weight * sqrtg;
                                                }
                                            }
                                        }
                                    }

                                }
                            }
                        }
                    }
                }
                break;

                // NOTE: Neumann Boundary conditions prescribe a solution gradient 
                // For this we only use the diffusive flux 
                // DiffusiveFlux must provide a neumann_flux function that takes the normal gradient
                // If the diffusive flux uses the solution value this will not work 
                // this also ignores the convective flux because hyperbolic problems 
                // dont have a notion of Neumann BC (use an outflow or extrapolation BC instead)
                case BOUNDARY_CONDITIONS::NEUMANN:
                {
                    // loop over quadrature points 
                    for(int iqp = 0; iqp < trace.nQP(); ++iqp){

                        const QuadraturePoint<T, ndim - 1> &quadpt = trace.getQP(iqp);

                        // calculate the jacobian and riemannian metric root det
                        auto Jfac = trace.face->Jacobian(coord, quadpt.abscisse);
                        T sqrtg = trace.face->rootRiemannMetric(Jfac, quadpt.abscisse);

                        // calculate the physical domain position
                        MATH::GEOMETRY::Point<T, ndim> phys_pt;
                        trace.face->transform(quadpt.abscisse, coord, phys_pt);

                        // get the basis function values
                        auto biL = trace.qp_evals_l[iqp].bi_span;

                        // Get the values at the boundary 
                        std::array<T, nv_comp> neumann_vals{};
                        neumann_callbacks[trace.face->bcflag](phys_pt.data(), neumann_vals.data());
                        // flux contribution weighted by quadrature and face metric 
                        // Li and Tang 2017 sec 9.1.1
                        std::array<T, nv_comp> fviscn = diff_flux.neumann_flux(neumann_vals);

                        // scale by weight and face metric tensor
                        for(int ieq = 0; ieq < neq; ++ieq){
                            fviscn[ieq] *= quadpt.weight * sqrtg;
                        }

                        // scatter contribution 
                        for(int itest = 0; itest < elL.nbasis(); ++itest){
                            for(int ieq = 0; ieq < neq; ++ieq){
                                resL[itest, ieq] += (fviscn[ieq]) * biL[itest];
                            }
                        }
                    }
                }
                break;

                case BOUNDARY_CONDITIONS::SPACETIME_PAST:
                if constexpr(!std::same_as<ST_Info, std::false_type>){

                    static constexpr int neq = nv_comp;
                    static_assert(neq == decltype(unkelL)::static_extent(), "Number of equations must match.");
                    static_assert(neq == decltype(unkelR)::static_extent(), "Number of equations must match.");
                    using namespace NUMTOOL::TENSOR::FIXED_SIZE;

                    // Get the info from the connection
                    TraceSpace<T, IDX, ndim>& trace_past = spacetime_info.connection_traces[trace.facidx];
                    std::vector<T> unker_past_data{trace_past.elL.nbasis() * neq};
                    auto uR_layout = spacetime_info.u_past.create_element_layout(trace_past.elR.elidx);
                    dofspan unkel_past{unker_past_data, uR_layout};
                    extract_elspan(trace_past.elR.elidx, spacetime_info.u_past, unkel_past);
                    AbstractMesh<T, IDX, ndim>* mesh_past = spacetime_info.fespace_past.meshptr;

                    // calculate the centroids of the left and right elements
                    // in the physical domain
                    const FiniteElement &elL = trace.elL;
                    const FiniteElement &elR = trace_past.elR;
                    auto centroidL = elL.geo_el->centroid(coord);

                    // Basis function scratch space 
                    PhysDomainEvalStorage storageL{elL};
                    PhysDomainEvalStorage storageR{elR};

                    // solution scratch space 
                    std::array<T, neq> uL;
                    std::array<T, neq> uR;
                    std::array<T, neq * ndim> graduL_data;
                    std::array<T, neq * ndim> graduR_data;
                    std::array<T, neq * ndim> grad_ddg_data;
                    std::array<T, neq * ndim * ndim> hessuL_data;
                    std::array<T, neq * ndim * ndim> hessuR_data;


                    // loop over the quadrature points 
                    for(int iqp = 0; iqp < trace.nQP(); ++iqp){
                        const QuadraturePoint<T, ndim - 1> &quadpt = trace.getQP(iqp);

                        // calculate the riemannian metric tensor root
                        auto Jfac = trace.face->Jacobian(coord, quadpt.abscisse);
                        T sqrtg = trace.face->rootRiemannMetric(Jfac, quadpt.abscisse);

                        // calculate the normal vector 
                        auto normal = calc_ortho(Jfac);
                        auto unit_normal = normalize(normal);

                        // get the basis functions, derivatives, and hessians
                        // (derivatives are wrt the physical domain)
                        auto biL = trace.qp_evals_l[iqp].bi_span;
                        auto biR = trace.qp_evals_r[iqp].bi_span;
                        auto xiL = trace.transform_xiL(quadpt.abscisse);
                        auto xiR = trace.transform_xiR(quadpt.abscisse);
                        PhysDomainEval evalL{storageL, elL, xiL, trace.qp_evals_l[iqp]};
                        PhysDomainEval evalR{storageR, elR, xiR, trace.qp_evals_r[iqp]};

                        // construct the solution on the left and right
                        std::ranges::fill(uL, 0.0);
                        std::ranges::fill(uR, 0.0);
                        for(int ieq = 0; ieq < neq; ++ieq){
                            for(int ibasis = 0; ibasis < elL.nbasis(); ++ibasis)
                                { uL[ieq] += unkelL[ibasis, ieq] * biL[ibasis]; }
                            for(int ibasis = 0; ibasis < elR.nbasis(); ++ibasis)
                                { uR[ieq] += unkel_past[ibasis, ieq] * biR[ibasis]; }
                        }

                        // get the solution gradient and hessians
                        auto graduL = unkelL.contract_mdspan(evalL.phys_grad_basis, graduL_data.data());
                        auto graduR = unkelR.contract_mdspan(evalR.phys_grad_basis, graduR_data.data());
                        auto hessuL = unkelL.contract_mdspan(evalL.phys_hess_basis, hessuL_data.data());
                        auto hessuR = unkelR.contract_mdspan(evalR.phys_hess_basis, hessuR_data.data());

                        // compute convective fluxes
                        std::array<T, neq> fadvn = conv_nflux(uL, uR, unit_normal);

                        // compute a single valued gradient using DDG or IP 

                        // calculate the DDG distance
                        MATH::GEOMETRY::Point<T, ndim> phys_pt;
                        trace.face->transform(quadpt.abscisse, coord, phys_pt);
                        T h_ddg = 0;
                        for(int idim = 0; idim < ndim; ++idim){
                            h_ddg += unit_normal[idim] * (
                                2 * (phys_pt[idim] - centroidL[idim])
                            );
                        }
                        h_ddg = std::copysign(std::max(std::abs(h_ddg), std::numeric_limits<T>::epsilon()), h_ddg);
                        
                        int order = std::max(
                            elL.basis->getPolynomialOrder(),
                            elR.basis->getPolynomialOrder()
                        );
                        // Danis and Yan reccomended for NS
                        T beta0 = std::pow(order + 1, 2);
                        T beta1 = 1 / std::max((T) (2 * order * (order + 1)), 1.0);

                        // switch to interior penalty if set
                        if(interior_penalty) beta1 = 0.0;


                        std::mdspan<T, std::extents<int, neq, ndim>> grad_ddg{grad_ddg_data.data()};
                        for(int ieq = 0; ieq < neq; ++ieq){
                            // construct the DDG derivatives
                            T jumpu = uR[ieq] - uL[ieq];
                            for(int idim = 0; idim < ndim; ++idim){
                                grad_ddg[ieq, idim] = beta0 * jumpu / h_ddg * unit_normal[idim]
                                    + 0.5 * (graduL[ieq, idim] + graduR[ieq, idim]);
                                T hessTerm = 0;
                                for(int jdim = 0; jdim < ndim; ++jdim){
                                    hessTerm += (hessuR[ieq, jdim, idim] - hessuL[ieq, jdim, idim])
                                        * unit_normal[jdim];
                                }
                                grad_ddg[ieq, idim] += beta1 * h_ddg * hessTerm;
                            }
                        }

                        // construct the viscous fluxes 
                        std::array<T, neq> uavg;
                        for(int ieq = 0; ieq < neq; ++ieq) uavg[ieq] = 0.5 * (uL[ieq] + uR[ieq]);

                        std::array<T, neq> fviscn = diff_flux(uavg, grad_ddg, unit_normal);

                        // scale by weight and face metric tensor
                        for(int ieq = 0; ieq < neq; ++ieq){
                            fadvn[ieq] *= quadpt.weight * sqrtg;
                            fviscn[ieq] *= quadpt.weight * sqrtg;
                        }

                        // scatter contribution 
                        for(int itest = 0; itest < elL.nbasis(); ++itest){
                            for(int ieq = 0; ieq < neq; ++ieq){
                                resL[itest, ieq] += (fviscn[ieq] - fadvn[ieq]) * biL[itest];
                            }
                        }
                    }
                    
                }
                break;

                case BOUNDARY_CONDITIONS::SPACETIME_FUTURE:
                    // NOTE: Purely upwind so continue on to use EXTRAPOLATION

                // Use only the interior state and assume the exterior state (and gradients) match
                case BOUNDARY_CONDITIONS::EXTRAPOLATION:
                {
                    // loop over quadrature points
                    for(int iqp = 0; iqp < trace.nQP(); ++iqp){
                        const QuadraturePoint<T, ndim - 1> &quadpt = trace.getQP(iqp);

                        // calculate the jacobian and riemannian metric root det
                        auto Jfac = trace.face->Jacobian(coord, quadpt.abscisse);
                        T sqrtg = trace.face->rootRiemannMetric(Jfac, quadpt.abscisse);

                        // calculate the normal vector 
                        auto normal = calc_ortho(Jfac);
                        auto unit_normal = normalize(normal);

                        // calculate the physical domain position
                        MATH::GEOMETRY::Point<T, ndim> phys_pt;
                        trace.face->transform(quadpt.abscisse, coord, phys_pt);

                        // get the function values
                        auto biL = trace.qp_evals_l[iqp].bi_span;

                        // get the gradients the physical domain
                        auto gradBiL = trace.eval_phys_grad_basis_l_qp(iqp, gradbL_data.data());

                        auto graduL = unkelL.contract_mdspan(gradBiL, graduL_data.data());

                        // construct the solution on the left and right
                        std::ranges::fill(uL, 0.0);
                        for(int ieq = 0; ieq < neq; ++ieq){
                            for(int ibasis = 0; ibasis < elL.nbasis(); ++ibasis)
                                { uL[ieq] += unkelL[ibasis, ieq] * biL[ibasis]; }
                        }

                        // compute convective fluxes
                        std::array<T, neq> fadvn = conv_nflux(uL, uL, unit_normal);

                        // construct the DDG derivatives
                        std::mdspan<T, std::extents<int, neq, ndim>> grad_ddg{grad_ddg_data.data()};
                        for(int ieq = 0; ieq < neq; ++ieq){
                            // construct the DDG derivatives ( just match interior gradient )
                            for(int idim = 0; idim < ndim; ++idim){
                                grad_ddg[ieq, idim] = (graduL[ieq, idim]);
                            }
                        }

                        // construct the viscous fluxes 
                        std::array<T, neq> uavg;
                        for(int ieq = 0; ieq < neq; ++ieq) uavg[ieq] = uL[ieq];

                        std::array<T, neq> fviscn = diff_flux(uavg, grad_ddg, unit_normal);

                        // scale by weight and face metric tensor
                        for(int ieq = 0; ieq < neq; ++ieq){
                            fadvn[ieq] *= quadpt.weight * sqrtg;
                            fviscn[ieq] *= quadpt.weight * sqrtg;
                        }

                        // scatter contribution 
                        for(int itest = 0; itest < elL.nbasis(); ++itest){
                            for(int ieq = 0; ieq < neq; ++ieq){
                                resL[itest, ieq] += (fviscn[ieq] - fadvn[ieq]) * biL[itest];
                            }
                        }
                    }
                }
                break;

                default:
                // === General BC Case ===
                // We construct the interior state uL 
                // and the gradients of the interior state graduL 
                // and pass this to the PDE implementation to give 
                // uR and graduR
                // otherwise
                // assume essential
                if constexpr (implements_bcs<PFlux>){
                    // loop over quadrature points 
                    for(int iqp = 0; iqp < trace.nQP(); ++iqp){
                        const QuadraturePoint<T, ndim - 1> &quadpt = trace.getQP(iqp);

                        // calculate the jacobian and riemannian metric root det
                        auto Jfac = trace.face->Jacobian(coord, quadpt.abscisse);
                        T sqrtg = trace.face->rootRiemannMetric(Jfac, quadpt.abscisse);

                        // calculate the normal vector 
                        auto normal = calc_ortho(Jfac);
                        auto unit_normal = normalize(normal);

                        // calculate the physical domain position
                        MATH::GEOMETRY::Point<T, ndim> phys_pt;
                        trace.face->transform(quadpt.abscisse, coord, phys_pt);

                        // get the function values
                        auto biL = trace.qp_evals_l[iqp].bi_span;

                        // construct the solution on the left
                        std::ranges::fill(uL, 0.0);
                        for(int ieq = 0; ieq < neq; ++ieq){
                            for(int ibasis = 0; ibasis < elL.nbasis(); ++ibasis)
                                { uL[ieq] += unkelL[ibasis, ieq] * biL[ibasis]; }
                        }

                        // get the gradients the physical domain
                        auto gradBiL = trace.eval_phys_grad_basis_l_qp(iqp, gradbL_data.data());
                        auto graduL = unkelL.contract_mdspan(gradBiL, graduL_data.data());

                        // get the uR and graduR from the interior information
                        // uL is also modifiable
                        auto [uR, graduR] = phys_flux.apply_bc(
                                uL, graduL, unit_normal,
                                trace.face->bctype, trace.face->bcflag);

                        // construct the DDG derivatives
                        int order = 
                            elL.basis->getPolynomialOrder();
                        // Danis and Yan reccomended for NS
                        T beta0 = std::pow(order + 1, 2);
                        T beta1 = 1 / std::max((T) (2 * order * (order + 1)), 1.0);

                        // switch to interior penalty if set
                        if(interior_penalty) beta1 = 0.0;

                        // calculate the DDG distance
                        T h_ddg = 0; // uses distance to quadpt on boundary face
                        for(int idim = 0; idim < ndim; ++idim){
                            h_ddg += std::abs(unit_normal[idim] * 
                                (phys_pt[idim] - centroidL[idim])
                            );
                        }
                        h_ddg = std::copysign(std::max(std::abs(h_ddg)
                                    , std::numeric_limits<T>::epsilon()), h_ddg);
                        std::mdspan<T, std::extents<int, neq, ndim>> 
                            grad_ddg{grad_ddg_data.data()};

                        for(int ieq = 0; ieq < neq; ++ieq){
                            // construct the DDG derivatives
                            T jumpu = uR[ieq] - uL[ieq];
                            for(int idim = 0; idim < ndim; ++idim){
                                grad_ddg[ieq, idim] = beta0 * jumpu / h_ddg * unit_normal[idim]
                                    + (graduL[ieq, idim]);
                            }
                        }

                        // compute convective fluxes
                        std::array<T, neq> fadvn = conv_nflux(uL, uR, unit_normal);

                        // construct the viscous fluxes 
                        std::array<T, neq> uavg;
                        for(int ieq = 0; ieq < neq; ++ieq) uavg[ieq] = 0.5 * (uL[ieq] + uR[ieq]);
                        std::array<T, neq> fviscn = diff_flux(uavg, grad_ddg, unit_normal);

                        // scale by weight and face metric tensor
                        for(int ieq = 0; ieq < neq; ++ieq){
                            fadvn[ieq] *= quadpt.weight * sqrtg;
                            fviscn[ieq] *= quadpt.weight * sqrtg;
                        }

                        // scatter contribution 
                        for(int itest = 0; itest < elL.nbasis(); ++itest){
                            for(int ieq = 0; ieq < neq; ++ieq){
                                resL[itest, ieq] += (fviscn[ieq] - fadvn[ieq]) * biL[itest];
                            }
                        }

                        // if applicable: apply the interface correction 
                        if constexpr (computes_homogeneity_tensor<DiffusiveFlux>) {
                            if(sigma_ic != 0.0){

                                std::array<T, neq> interface_correction;
                                auto Gtensor = diff_flux.homogeneity_tensor(uavg);

                                T average_gradv[ndim];
                                for(int itest = 0; itest < elL.nbasis(); ++itest){
                                    // get the average test function gradient
                                    for(int idim = 0; idim < ndim; ++idim){
                                        average_gradv[idim] = gradBiL[itest, idim];
                                    }

                                    std::ranges::fill(interface_correction, 0);
                                    for(int ieq = 0; ieq < neq; ++ieq){
                                        for(int kdim = 0; kdim < ndim; ++kdim){
                                            for(int req = 0; req < neq; ++req){
                                                T jumpu_r = uR[req] - uL[req];
                                                for(int sdim = 0; sdim < ndim; ++sdim){
                                                    resL[itest, ieq] -= 
                                                        Gtensor[ieq][kdim][req][sdim] * unit_normal[kdim] 
                                                        * average_gradv[sdim] * jumpu_r
                                                        * quadpt.weight * sqrtg;
                                                }
                                            }
                                        }
                                    }

                                }
                            }
                        }
                    }
                }
                break;
            }
        }

        template<class IDX>
        void interface_conservation(
            const TraceSpace<T, IDX, ndim>& trace,
            NodeArray<T, ndim>& coord,
            elspan auto unkelL,
            elspan auto unkelR,
            facspan auto res
        ) const {
            static constexpr int neq = nv_comp;
            using namespace MATH::MATRIX_T;
            using namespace NUMTOOL::TENSOR::FIXED_SIZE;
            using FiniteElement = FiniteElement<T, IDX, ndim>;

            const FiniteElement &elL = trace.elL;
            const FiniteElement &elR = trace.elR;

            // Basis function scratch space 
            std::vector<T> bitrace(trace.nbasis_trace());
            std::vector<T> gradbL_data(elL.nbasis() * ndim);
            std::vector<T> gradbR_data(elR.nbasis() * ndim);

            // solution scratch space 
            std::array<T, neq> uL;
            std::array<T, neq> uR;
            std::array<T, neq * ndim> graduL_data;
            std::array<T, neq * ndim> graduR_data;
            std::array<T, neq * ndim> grad_ddg_data;

            for(int iqp = 0; iqp < trace.nQP(); ++iqp){
                const QuadraturePoint<T, ndim - 1> &quadpt = trace.getQP(iqp);

                // calculate the riemannian metric tensor root
                auto Jfac = trace.face->Jacobian(coord, quadpt.abscisse);
                T sqrtg = trace.face->rootRiemannMetric(Jfac, quadpt.abscisse);

                // calculate the normal vector 
                auto normal = calc_ortho(Jfac);
                auto unit_normal = normalize(normal);

                // get the basis functions, derivatives, and hessians
                // (derivatives are wrt the physical domain)
                auto biL = trace.eval_basis_l_qp(iqp);
                auto biR = trace.eval_basis_r_qp(iqp);
                trace.eval_trace_basis_qp(iqp, bitrace.data());
                auto gradBiL = trace.eval_phys_grad_basis_l_qp(iqp, gradbL_data.data());
                auto gradBiR = trace.eval_phys_grad_basis_r_qp(iqp, gradbR_data.data());

                // construct the solution on the left and right
                std::ranges::fill(uL, 0.0);
                std::ranges::fill(uR, 0.0);
                for(int ieq = 0; ieq < neq; ++ieq){
                    for(int ibasis = 0; ibasis < elL.nbasis(); ++ibasis)
                        { uL[ieq] += unkelL[ibasis, ieq] * biL[ibasis]; }
                    for(int ibasis = 0; ibasis < elR.nbasis(); ++ibasis)
                        { uR[ieq] += unkelR[ibasis, ieq] * biR[ibasis]; }
                }

                // get the solution gradient and hessians
                auto graduL = unkelL.contract_mdspan(gradBiL, graduL_data.data());
                auto graduR = unkelR.contract_mdspan(gradBiR, graduR_data.data());

                if(trace.face->bctype != BOUNDARY_CONDITIONS::INTERIOR){
                    switch(trace.face->bctype){
                        case BOUNDARY_CONDITIONS::DIRICHLET:
                            {
                                if(trace.elL.elidx != trace.elR.elidx)
                                    std::cout << "warning: elements do not match" << std::endl;

                                // calculate the physical domain position
                                MATH::GEOMETRY::Point<T, ndim> phys_pt;
                                trace.face->transform(quadpt.abscisse, coord, phys_pt);

                                // Get the values at the boundary 
                                dirichlet_callbacks[trace.face->bcflag](phys_pt.data(), uR.data());

                            } break;
                        default:
                            for(int itest = 0; itest < trace.nbasis_trace(); ++itest){
                                for(int ieq = 0; ieq < neq; ++ieq)
                                    res[itest, ieq] = 0;
                            } return;
                    }
                }

                // HACK: disable diffusion IC for linear polynomials 
                if(elL.basis->getPolynomialOrder() == 1 && elR.basis->getPolynomialOrder() == 1){
                    std::ranges::fill(graduL_data, 0.0);
                    std::ranges::fill(graduR_data, 0.0);
                }

                // get the physical flux on the left and right
                Tensor<T, neq, ndim> fluxL = phys_flux(uL, graduL);
                Tensor<T, neq, ndim> fluxR = phys_flux(uR, graduR);

                // calculate the jump in normal fluxes
                Tensor<T, neq> jumpflux{};
                for(int ieq = 0; ieq < neq; ++ieq) 
                    jumpflux[ieq] = dot(fluxR[ieq], unit_normal) - dot(fluxL[ieq], unit_normal);

                // scatter unit normal times interface conservation to residual
                for(int itest = 0; itest < trace.nbasis_trace(); ++itest){
                    for(int ieq = 0; ieq < neq; ++ieq){
                        T ic_res = jumpflux[ieq] * sqrtg * quadpt.weight;
                        // NOTE: multiplying by signed unit normal 
                        // adds directionality which can allow cancellation error with 
                        // V-shaped interface intersections
                        res[itest, ieq] -= ic_res * bitrace[itest]; 
                    }
                }
            }
        }
    };

    // Deduction Guides
    template<
        class PFlux,
        class CFlux, 
        class DiffusiveFlux
    >
    ConservationLawDDG(
        PFlux&& physical_flux,
        CFlux&& convective_numflux,
        DiffusiveFlux&& diffusive_flux
    ) -> ConservationLawDDG<typename PFlux::value_type, PFlux::ndim, PFlux, CFlux, DiffusiveFlux>;

    template<
        class PFlux,
        class CFlux, 
        class DiffusiveFlux,
        class ST_Info
    >
    ConservationLawDDG(
        PFlux&& physical_flux,
        CFlux&& convective_numflux,
        DiffusiveFlux&& diffusive_flux,
        ST_Info st_info
    ) -> ConservationLawDDG<typename PFlux::value_type, PFlux::ndim, PFlux, CFlux, DiffusiveFlux, ST_Info>;
}
