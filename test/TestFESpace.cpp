#include "iceicle/disc/projection.hpp"
#include "iceicle/element/finite_element.hpp"
#include "iceicle/element/reference_element.hpp"
#include "iceicle/fe_enums.hpp"
#include "iceicle/fe_function/dglayout.hpp"
#include "iceicle/fe_function/el_layout.hpp"
#include "iceicle/geometry/geo_element.hpp"
#include "iceicle/mesh/mesh.hpp"
#include "iceicle/quadrature/HypercubeGaussLegendre.hpp"
#include "iceicle/quadrature/quadrules_1d.hpp"
#include "iceicle/solvers/element_linear_solve.hpp"
#include "iceicle/tmp_utils.hpp"
#include <iceicle/fespace/fespace.hpp>
#include <iceicle/fe_utils.hpp>

#include <gtest/gtest.h>
#include <pstl/glue_execution_defs.h>

TEST(test_fespace, test_element_construction){
    using T = double;
    using IDX = int;

    static constexpr int ndim = 2;
    static constexpr int pn_geo = 2;
    static constexpr int pn_basis = 3;

    // create a uniform mesh
    MESH::AbstractMesh<T, IDX, ndim> mesh({-1.0, -1.0}, {1.0, 1.0}, {2, 2}, pn_basis);

    FE::FESpace<T, IDX, ndim> fespace{
        &mesh, FE::FESPACE_ENUMS::LAGRANGE,
        FE::FESPACE_ENUMS::GAUSS_LEGENDRE, 
        ICEICLE::TMP::compile_int<pn_basis>()
    };

    ASSERT_EQ(fespace.elements.size(), 4);

    ASSERT_EQ(fespace.dg_offsets.calculate_size_requirement(2), 4 * 2 * std::pow(pn_basis + 1, pn_geo));
}

class test_geo_el : public ELEMENT::GeometricElement<double, int, 2>{
    static constexpr int ndim = 2;
    using T = double;
    using IDX = int;

    using Point = MATH::GEOMETRY::Point<T, ndim>;
    using HessianType = NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim, ndim>;

public:
    
    constexpr int n_nodes() const override {return 0;};

    constexpr FE::DOMAIN_TYPE domain_type() const noexcept override {return FE::DOMAIN_TYPE::DYNAMIC;};

    constexpr int geometry_order() const noexcept override {return 1;};

    const IDX *nodes() const override { return nullptr; }

    void transform(
        FE::NodalFEFunction<T, ndim> &node_coords,
        const Point &pt_ref,
        Point &pt_phys
    ) const override {
        T xi = pt_ref[0];
        T eta = pt_ref[1];

        pt_phys[0] = xi * eta;
        pt_phys[1] = xi + eta;
    }

    NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim> Jacobian(
        FE::NodalFEFunction<T, ndim> &node_coords,
        const Point &xi_arg
    ) const override {
        T xi = xi_arg[0];
        T eta = xi_arg[1];
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim> jac = {{
            {eta, xi},
            {1, 1}
        }};
        return jac;
    }

    HessianType Hessian(
        FE::NodalFEFunction<T, ndim> &node_coords,
        const Point &xi
    ) const override {
        HessianType hess;
        hess[0][0][0] = 0;
        hess[0][0][1] = 1;
        hess[0][1][0] = 1;
        hess[0][1][1] = 0;
        hess[1][0][0] = 0;
        hess[1][0][1] = 0;
        hess[1][1][0] = 0;
        hess[1][1][1] = 0;
        return hess;
    }
};

class test_basis : public BASIS::Basis<double, 2> {
    static constexpr int ndim = 2;
    using T = double;
public:

    int nbasis() const override {return 1;}

    constexpr FE::DOMAIN_TYPE domain_type() const noexcept override { return FE::DOMAIN_TYPE::DYNAMIC;}

    void evalBasis(const T *xi_vec, T *Bi) const override {
        T xi = xi_vec[0];
        T eta = xi_vec[1];

        Bi[0] = eta * xi * xi + eta * eta * xi;
    }

    void evalGradBasis(const T *xi_vec, T *dBidxj) const override {
        T xi = xi_vec[0];
        T eta = xi_vec[1];
        dBidxj[0] = 2 * xi * eta + eta * eta;
        dBidxj[1] = xi * xi + 2 * xi * eta;

    }

    void evalHessBasis(const T*xi_vec, T *Hessian) const override {
        T xi = xi_vec[0];
        T eta = xi_vec[1];
       
        Hessian[0] = 2 * eta;
        Hessian[1] = 2 * xi + 2 * eta;
        Hessian[2] = Hessian[1];
        Hessian[3] = 2 * xi;
    };

    bool isOrthonormal() const override { return false; }

    bool isNodal() const override { return false; }

    inline int getPolynomialOrder() const override {return 2;}

};

TEST(test_finite_element, test_hess_basis){
    static constexpr int ndim = 2;

    QUADRATURE::HypercubeGaussLegendre<double, int, 2, 1> unused{};
    test_basis basis{};
    test_geo_el geo_el{};

    FE::NodalFEFunction<double, 2> coord{};

    ELEMENT::FEEvaluation<double, int, 2> evals_unused{};

    ELEMENT::FiniteElement<double, int, 2> el{&geo_el, &basis, &unused, &evals_unused, 0};

    double xi = -0.2;
    double eta = 0.5;
    MATH::GEOMETRY::Point<double, 2> ref_pt = {xi, eta};

    // test he hessian
    std::vector<double> hess_basis_data(el.nbasis() * ndim * ndim);
    auto hess_basis = el.evalPhysHessBasis(ref_pt, coord, hess_basis_data.data());

    ASSERT_NEAR((hess_basis[0, 0, 0]), 0, 1e-14);
    ASSERT_NEAR((hess_basis[0, 0, 1]), 1, 1e-14);
    ASSERT_NEAR((hess_basis[0, 1, 0]), 1, 1e-14);
    ASSERT_NEAR((hess_basis[0, 1, 1]), 0, 1e-14);
}

TEST(test_fespace, test_dg_projection){

    using T = double;
    using IDX = int;

    static constexpr int ndim = 2;
    static constexpr int pn_geo = 1;
    static constexpr int pn_basis = 4;
    static constexpr int neq = 1;

    // create a uniform mesh
    int nx = 50;
    int ny = 10;
    MESH::AbstractMesh<T, IDX, ndim> mesh({-1.0, -1.0}, {1.0, 1.0}, {nx, ny}, pn_geo);
    mesh.nodes.random_perturb(-0.4 * 1.0 / std::max(nx, ny), 0.4*1.0/std::max(nx, ny));

    FE::FESpace<T, IDX, ndim> fespace{
        &mesh, FE::FESPACE_ENUMS::LAGRANGE,
        FE::FESPACE_ENUMS::GAUSS_LEGENDRE, 
        ICEICLE::TMP::compile_int<pn_basis>()
    };

    // define a pn_basis order polynomial function to project onto the space 
    auto projfunc = [](const double *xarr, double *out){
        double x = xarr[0];
        double y = xarr[1];
        out[0] = std::pow(x, pn_basis) + std::pow(y, pn_basis);
    };

    auto dprojfunc = [](const double *xarr) {
        double x = xarr[0];
        double y = xarr[1];
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<double, ndim> deriv = { 
            pn_basis * std::pow(x, pn_basis - 1),
            pn_basis * std::pow(y, pn_basis - 1)
        };
        return deriv;
    };

    auto hessfunc = [](const double *xarr){
        double x = xarr[0];
        double y = xarr[1];
        int n = pn_basis;
        if (pn_basis < 2){
            NUMTOOL::TENSOR::FIXED_SIZE::Tensor<double, ndim, ndim> hess = {{
                {0, 0},
                {0, 0}
            }};
            return hess;
        } else {
            NUMTOOL::TENSOR::FIXED_SIZE::Tensor<double, ndim, ndim> hess = {{
                {n * (n-1) *std::pow(x, n-2), 0.0},
                {0.0, n * (n-1) *std::pow(y, n-2)}
            }};
            return hess;
        }
    };


    // create the projection discretization
    DISC::Projection<double, int, ndim, neq> projection{projfunc};

    T *u = new T[fespace.ndof_dg() * neq](); // 0 initialized
    FE::fespan<T, FE::dg_layout<T, 1>> u_span(u, fespace.dg_offsets);

    // solve the projection 
    std::for_each(fespace.elements.begin(), fespace.elements.end(),
        [&](const ELEMENT::FiniteElement<T, IDX, ndim> &el){
            FE::compact_layout<double, 1> el_layout{el};
            T *u_local = new T[el_layout.size()](); // 0 initialized 
            FE::elspan<T, FE::compact_layout<double, 1>> u_local_span(u_local, el_layout);

            T *res_local = new T[el_layout.size()](); // 0 initialized 
            FE::elspan<T, FE::compact_layout<double, 1>> res_local_span(res_local, el_layout);
            
            // projection residual
            projection.domainIntegral(el, fespace.meshptr->nodes, res_local_span);

            // solve 
            SOLVERS::ElementLinearSolver<T, IDX, ndim, neq> solver{el, fespace.meshptr->nodes};
            solver.solve(u_local_span, res_local_span);

            // test a bunch of random locations
            for(int k = 0; k < 50; ++k){
                MATH::GEOMETRY::Point<T, ndim> ref_pt = FE::random_domain_point(el.geo_el);
                MATH::GEOMETRY::Point<T, ndim> phys_pt;
                el.transform(fespace.meshptr->nodes, ref_pt, phys_pt);
                
                // get the actual value of the function at the given point in the physical domain
                T act_val;
                projfunc(phys_pt, &act_val);

                T projected_val = 0;
                T *basis_vals = new double[el.nbasis()];
                el.evalBasis(ref_pt, basis_vals);
                u_local_span.contract_dofs(basis_vals, &projected_val);

                ASSERT_NEAR(projected_val, act_val, 1e-8);

                // test the derivatives
                std::vector<double> grad_basis_data(el.nbasis() * ndim);
                auto grad_basis = el.evalPhysGradBasis(ref_pt, fespace.meshptr->nodes, grad_basis_data.data());
                static_assert(grad_basis.rank() == 2);
                static_assert(grad_basis.extent(1) == ndim);

                // get the derivatives for each equation by contraction
                std::vector<double> grad_eq_data(neq * ndim, 0);
                auto grad_eq = u_local_span.contract_mdspan(grad_basis, grad_eq_data.data());

                auto dproj = dprojfunc(phys_pt);
                ASSERT_NEAR(dproj[0], (grad_eq[0, 0]), 1e-8);
                ASSERT_NEAR(dproj[1], (grad_eq[0, 1]), 1e-8);

                // test hessian
                std::vector<double> hess_basis_data(el.nbasis() * ndim * ndim);
                auto hess_basis = el.evalPhysHessBasis(ref_pt, fespace.meshptr->nodes, hess_basis_data.data());

                // get the hessian for each equation by contraction 
                std::vector<double> hess_eq_data(neq * ndim * ndim, 0);
                auto hess_eq = u_local_span.contract_mdspan(hess_basis, hess_eq_data.data());
                auto hess_proj = hessfunc(phys_pt);
                ASSERT_NEAR(hess_proj[0][0], (hess_eq[0, 0, 0]), 1e-8);
                ASSERT_NEAR(hess_proj[0][1], (hess_eq[0, 0, 1]), 1e-8);
                ASSERT_NEAR(hess_proj[1][0], (hess_eq[0, 1, 0]), 1e-8);
                ASSERT_NEAR(hess_proj[1][1], (hess_eq[0, 1, 1]), 1e-8);

                delete[] basis_vals;
            }

            delete[] u_local;
            delete[] res_local;
        }
    );

    delete[] u;



}
