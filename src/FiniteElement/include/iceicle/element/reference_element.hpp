/**
 * @brief reference element stores things that are common between FiniteElements
 * @author Gianni Absillis (gabsill@ncsu.edu)
 */
#pragma once
#include "iceicle/basis/lagrange.hpp"
#include "iceicle/fe_enums.hpp"
#include "iceicle/geometry/geometry_enums.hpp"
#include "iceicle/quadrature/SimplexQuadrature.hpp"
#include "iceicle/quadrature/HypercubeGaussLegendre.hpp"
#include "iceicle/quadrature/quadrules_1d.hpp"
#include <iceicle/element/finite_element.hpp>
#include <iceicle/tmp_utils.hpp>
#include <memory>

namespace FE {
    namespace FESPACE_ENUMS {

        enum FESPACE_BASIS_TYPE {
            /// Lagrange polynomials
            LAGRANGE = 0, 
            N_BASIS_TYPES,
        };

        enum FESPACE_QUADRATURE {
            /// Gauss Legendre Quadrature rules and tensor extensions thereof
            /// uses Grundmann Moller for Simplex type elements
            GAUSS_LEGENDRE,
            N_QUADRATURE_TYPES
        };
    }
}

namespace ELEMENT {

    template<typename T, typename IDX, int ndim>
    class ReferenceElement {
        using FEEvalType = ELEMENT::FEEvaluation<T, IDX, ndim>;
        using BasisType = BASIS::Basis<T, ndim>;
        using QuadratureType = QUADRATURE::QuadratureRule<T, IDX, ndim>;
        using GeoElementType = ELEMENT::GeometricElement<T, IDX, ndim>;

        public:
        std::unique_ptr<BasisType> basis;
        std::unique_ptr<QuadratureType> quadrule;
        FEEvalType eval;

        ReferenceElement() = default;

        template<int basis_order>
        ReferenceElement(
            const GeoElementType *geo_el,
            FE::FESPACE_ENUMS::FESPACE_BASIS_TYPE basis_type,
            FE::FESPACE_ENUMS::FESPACE_QUADRATURE quadrature_type,
            ICEICLE::TMP::compile_int<basis_order> basis_order_arg
        ) {
            using namespace FE;
            using namespace FESPACE_ENUMS;
            switch(geo_el->domain_type()){
                case FE::HYPERCUBE:
                    // construct the basis 
                    switch(basis_type){
                        case LAGRANGE:
                            basis = std::make_unique<BASIS::HypercubeLagrangeBasis<
                                T, IDX, ndim, basis_order>>();
                            break;
                        default:
                            break;
                    }

                    // construct the quadrature rule
                    // TODO: change quadrature order based on high order geo elements
                    switch(quadrature_type){
                        case FE::FESPACE_ENUMS::GAUSS_LEGENDRE:
                            quadrule = std::make_unique<QUADRATURE::HypercubeGaussLegendre<T, IDX, ndim, basis_order+1>>();
                            break;
                        default:
                            break;
                    }

                    // construct the evaluation
                    eval = FEEvalType(basis.get(), quadrule.get());

                    break;

                case FE::SIMPLEX:
                    // construct the basis 
                    switch(basis_type){
                        case LAGRANGE:
                            basis = std::make_unique<BASIS::SimplexLagrangeBasis<
                                T, IDX, ndim, basis_order>>();
                            break;
                        default:
                            break;
                    }

                    // construct the quadrature rule
                    switch(quadrature_type){
                        case FE::FESPACE_ENUMS::GAUSS_LEGENDRE:
                            quadrule = std::make_unique<QUADRATURE::GrundmannMollerSimplexQuadrature<T, IDX, ndim, basis_order+1>>();
                            break;
                        default:
                            break;
                    }

                    // construct the evaluation
                    eval = FEEvalType(basis.get(), quadrule.get());

                    break;
                default: 
                    break;
            }
        }
        
    };
}