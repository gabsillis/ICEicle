/// @brief Levenberg-Marquardt implementation following Ching et al. 
/// The moving discontinuous Galerkin method with interface condition enforcement for the simulation of hypersonic, viscous flows
/// Computer Methods in Applied Mechanics and Engineering 2024
///
/// @author Gianni Absillis (gabsill@ncsu.edu)

#include "Numtool/fixed_size_tensor.hpp"
#include "iceicle/anomaly_log.hpp"
#include "iceicle/fe_function/component_span.hpp"
#include "iceicle/fe_function/fespan.hpp"
#include "iceicle/fe_function/geo_layouts.hpp"
#include "iceicle/fe_function/layout_right.hpp"
#include "iceicle/fespace/fespace.hpp"
#include "iceicle/form_dense_jacobian.hpp"
#include "iceicle/nonlinear_solver_utils.hpp"
#include "iceicle/petsc_interface.hpp"
#include "iceicle/form_petsc_jacobian.hpp"
#include "iceicle/form_residual.hpp"

#include <petsc.h>

#include <iostream>
#include <petscmat.h>
#include <petscsys.h>
#include <petscvec.h>

namespace iceicle::solvers {

    namespace impl {
        /// @brief Petsc context for the Gauss Newton subproblem
        struct GNSubproblemCtx {
            /// the Jacobian Matrix
            Mat J; 

            /// The vector to use as temporary storage for Jx
            ///
            /// WARNING: must be setup by user
            Vec Jx;

            /// @brief regularization for pde dofs
            PetscScalar lambda_u = 1e-7;

            /// @brief anisotropic lagrangian regularization 
            PetscScalar lambda_lag = 1e-5;

            /// @brief Curvature penalization
            PetscScalar lambda_1 = 1e-3;

            /// @brief Grid pentalty regularization 
            PetscScalar lambda_b = 1e-2;

            /// @brief power of anisotropic metric
            PetscScalar alpha = -1;

            /// @brief power for principle stretching magnitude
            PetscScalar beta = 3;

            /// @brief the minimum allowable jacobian determinant
            PetscScalar J_min = 1e-10;

            PetscInt npde;
            PetscInt ngeo;
        };


        inline
        auto gn_subproblem(Mat A, Vec x, Vec y) -> PetscErrorCode {
            GNSubproblemCtx *ctx;

            PetscFunctionBeginUser;
            // Allocate the storage for x1

            PetscCall(MatShellGetContext(A, &ctx));

            // Jx = J*x
            PetscCall(MatMult(ctx->J, x, ctx->Jx));

            // y = J^T*J*x
            PetscCall(MatMultTranspose(ctx->J, ctx->Jx, y));

            Vec lambdax;
            VecCreate(PETSC_COMM_WORLD, &lambdax);
            VecSetSizes(lambdax, ctx->npde + ctx->ngeo, PETSC_DETERMINE);
            VecSetFromOptions(lambdax);
            PetscScalar* lambdax_data;
            const PetscScalar* xdata;
            VecGetArray(lambdax, &lambdax_data);
            VecGetArrayRead(x, &xdata);

//          // TODO: switch to parallel full size
            std::vector<PetscScalar> colnorms(ctx->npde + ctx->ngeo);
            // scaling using column norms (More 1977 Levenberg-Marquardt Implementation and Theory)
            MatGetColumnNorms(ctx->J, NORM_2, colnorms.data());
            for(PetscInt i = 0; i < ctx->npde; ++i) {
                lambdax_data[i] = xdata[i] * ctx->lambda_u * colnorms[i];
            }
            for(PetscInt i = ctx->npde; i < ctx->npde + ctx->ngeo; ++i )
                { lambdax_data[i] = xdata[i] * std::max(ctx->lambda_b, ctx->lambda_b * colnorms[i]); }
            VecRestoreArray(lambdax, &lambdax_data);
            VecRestoreArrayRead(x, &xdata);

            // y = (J^T*J + lambda * I)*x
            PetscCall(VecAXPY(y, 1.0, lambdax));

            VecDestroy(&lambdax);
            PetscFunctionReturn(EXIT_SUCCESS);
        }
    }
    
    template<class T, class IDX, int ndim, class disc_class, class ls_type = no_linesearch<T, IDX>>
    class CorriganLM {
        public:

        
        // ================
        // = Data Members =
        // ================

        /// @brief reference to the fespace to use
        FESpace<T, IDX, ndim>& fespace;

        /// @brief reference to the discretization to use
        disc_class& disc;

        /// @brief the convergence crieria
        ///
        /// determines whether the solver should terminate
        ConvergenceCriteria<T, IDX>& conv_criteria;

        /// @brief the linesearch strategy
        const ls_type& linesearch;

        /// @brief map of geometry dofs to consider for interface conservation enforcement
        const geo_dof_map<T, IDX, ndim>& geo_map;

        // === Petsc Data Members ===

        /// @brief the Jacobian Matrix
        Mat jac;

        /// @brief the matrix for the linear subproblem (JTJ + Regularization)
        Mat subproblem_mat;

        /// @brief context for matrix free subproblem implementation
        impl::GNSubproblemCtx subproblem_ctx;

        /// @brief the residual vector 
        Vec res_data;

        /// @brief the solution update
        Vec du_data;

        /// @brief J transpose times r 
        Vec Jtr;

        /// @brief J times x (for matrix free)
        Vec Jx;

        /// @brief krylov solver 
        KSP ksp;

        /// @brief preconditioner
        PC pc;

        // === Nonlinear Solver Behavior ===

        IDX verbosity = 0;

        /// @brief if this is a positive integer 
        /// Then the diagnostics callback will be called every ivis timesteps
        /// (k % ivis == 0)        
        IDX ivis = -1;

        /// @brief the callback function for visualization during solve()
        ///
        /// is given a reference to this when called 
        /// default is to print out a l2 norm of the residual data array
        /// Passes a reference to this, the current iteration number, the residual vector, and the du vector
        std::function<void(IDX, Vec, Vec)> vis_callback = []
            (IDX k, Vec res_data, Vec du_data)
        {
            T res_norm;
            PetscCallAbort(PETSC_COMM_WORLD, VecNorm(res_data, NORM_2, &res_norm));
            std::cout << std::format("itime: {:6d} | residual l1: {:16.8f}", k, res_norm) << std::endl;
        };

        /// @brief if this is a positive integer 
        /// Then the diagnostics callback will be called every idiag timesteps
        /// (k % idiag == 0)
        IDX idiag = -1;

        /// @brief diagnostics function 
        ///
        /// very minimal by default other options are defined in this header
        /// or a custom function can be made 
        ///
        /// Passes the current iteration number, the residual vector, and the du vector
        std::function<void(IDX, Vec, Vec)> diag_callback = []
            (IDX k, Vec res_data, Vec du_data)
        {
            int iproc;
            MPI_Comm_rank(PETSC_COMM_WORLD, &iproc); if(iproc == 0){
                std::cout << "Diagnostics for iteration: " << k << std::endl;
            }
            if(iproc == 0) std::cout << "Residual: " << std::endl;
            PetscCallAbort(PETSC_COMM_WORLD, VecView(res_data, PETSC_VIEWER_STDOUT_WORLD));
            if(iproc == 0) std::cout << std::endl << "du: " << std::endl;
            PetscCallAbort(PETSC_COMM_WORLD, VecView(du_data, PETSC_VIEWER_STDOUT_WORLD));
            if(iproc == 0) std::cout << "------------------------------------------" << std::endl << std::endl; 
        };

        /// @brief set to true if you want to explicitly form the J^TJ + lambda * I matrix
        ///
        /// This may greatly reduce the sparsity
        const bool explicitly_form_subproblem;

        // === Regularization Parameters ===
        
        /// @brief regularization for pde dofs
        T lambda_u = 1e-7;

        /// @brief anisotropic lagrangian regularization 
        T lambda_lag = 1e-5;

        /// @brief Curvature penalization
        T lambda_1 = 1e-3;

        /// @brief Grid pentalty regularization 
        T lambda_b = 1e-2;

        /// @brief power of anisotropic metric
        T alpha = -1;

        /// @brief power for principle stretching magnitude
        T beta = 3;

        /// @brief the minimum allowable jacobian determinant
        T J_min = 1e-10;

        public:

        // ================
        // = Constructors =
        // ================
        
        CorriganLM(
            FESpace<T, IDX, ndim>& fespace,
            disc_class& disc, 
            ConvergenceCriteria<T, IDX>& conv_criteria,
            const ls_type& linesearch,
            const geo_dof_map<T, IDX, ndim>& geo_map,
            bool explicitly_form_subproblem = false
        ) : fespace{fespace}, disc{disc}, conv_criteria{conv_criteria}, 
            linesearch{linesearch}, geo_map{geo_map},
            explicitly_form_subproblem{explicitly_form_subproblem}
        {
            static constexpr int neq = disc_class::nv_comp;

            // define data layouts
            fe_layout_right u_layout{fespace.dg_map, std::integral_constant<std::size_t, neq>{}};
            geo_data_layout geo_layout{geo_map};
            ic_residual_layout<T, IDX, ndim, neq> ic_layout{geo_map};

            // determine the system sizes on the local processor
            PetscInt local_u_size = u_layout.size() + geo_layout.size();
            PetscInt local_res_size = u_layout.size() + ic_layout.size();

            // Create and set up the jacobian matrix 
            MatCreate(PETSC_COMM_WORLD, &jac);
            MatSetSizes(jac, local_res_size, local_u_size, PETSC_DETERMINE, PETSC_DETERMINE);
            MatSetFromOptions(jac);
            MatSetUp(jac);

            if(explicitly_form_subproblem){

                // create the product matrix
                MatProductCreate(jac, jac, NULL, &subproblem_mat);
                MatProductSetType(subproblem_mat, MATPRODUCT_AtB);
                MatProductSetFromOptions(subproblem_mat);

            } else {
                
                MatCreate(PETSC_COMM_WORLD, &subproblem_mat);
                MatSetSizes(subproblem_mat, local_u_size, local_u_size, PETSC_DETERMINE, PETSC_DETERMINE);
                // setup the context and Shell matrix
                VecCreate(PETSC_COMM_WORLD, &Jx);
                VecSetSizes(Jx, local_res_size, PETSC_DETERMINE);
                VecSetFromOptions(Jx);

                subproblem_ctx.J = jac;
                subproblem_ctx.Jx = Jx;
                subproblem_ctx.npde = fespace.dg_map.calculate_size_requirement(disc_class::nv_comp);
                subproblem_ctx.ngeo = geo_map.size();

                MatSetType(subproblem_mat, MATSHELL);
                MatSetUp(subproblem_mat);
                MatShellSetOperation(subproblem_mat, MATOP_MULT, (void (*)()) impl::gn_subproblem);
                MatShellSetContext(subproblem_mat, (void *) &subproblem_ctx);
            }

            // Create and set up vectors
            // Create and set up the vectors
            VecCreate(PETSC_COMM_WORLD, &res_data);
            VecSetSizes(res_data, local_res_size, PETSC_DETERMINE);
            VecSetFromOptions(res_data);
            

            VecCreate(PETSC_COMM_WORLD, &Jtr);
            VecSetSizes(Jtr, local_u_size, PETSC_DETERMINE);
            VecSetFromOptions(Jtr);

            VecCreate(PETSC_COMM_WORLD, &du_data);
            VecSetSizes(du_data, local_u_size, PETSC_DETERMINE);
            VecSetFromOptions(du_data);

            // Create the linear solver and preconditioner
            PetscCallAbort(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &ksp));
            PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(ksp));

            // default preconditioner
            PetscCallAbort(PETSC_COMM_WORLD, KSPGetPC(ksp, &pc));

            if(explicitly_form_subproblem){
                PCSetType(pc, PCILU);
            } else {
                PCSetType(pc, PCNONE);
            }

            // Get user input (can override defaults set above)
            PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(ksp));
        }

        // ====================
        // = Member Functions =
        // ====================

        template<class uLayoutPolicy>
        auto solve(fespan<T, uLayoutPolicy> u) -> IDX {

            // update context from current state to make sure everything is synced
            subproblem_ctx.lambda_u = lambda_u;
            subproblem_ctx.lambda_b = lambda_b;

            static constexpr int neq = disc_class::nv_comp;

            // define data layouts
            fe_layout_right u_layout{fespace.dg_map, std::integral_constant<std::size_t, neq>{}};
            geo_data_layout geo_layout{geo_map};
            ic_residual_layout<T, IDX, ndim, neq> ic_layout{geo_map};

            // get the current coordinate data
            std::vector<T> coord_data(geo_layout.size());
            component_span coord{coord_data, geo_layout};
            extract_geospan(*(fespace.meshptr), coord);

            // get initial residual and jacobian 
            {
                petsc::VecSpan res_view{res_data};
                fespan res{res_view.data(), u.get_layout()};
                form_petsc_jacobian_fd(fespace, disc, u, res, jac);
                dofspan mdg_res{res_view.data() + u_layout.size(), ic_layout};
                form_petsc_mdg_jacobian_fd(fespace, disc, u, coord, mdg_res, jac);

                // assemble the Jacobian matrix  (assembly needed for symbolic product)
                MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
//
//                // Test Jacobian correctness
//                std::vector<T> densejac_data( (res.size() + mdg_res.size()) * (u.size() + coord.size()) );
//                std::mdspan<T, std::extents<std::size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_right> densejac{
//                    densejac_data.data(), res.size() + mdg_res.size(), u.size() + coord.size()
//                };
//                std::vector<T> testres(res.size() + mdg_res.size());
//                std::vector<T> testu(u.size() + coord.size());
//                std::copy_n(u.data(), u.size(), testu.begin());
//                std::copy_n(coord.data(), coord.size(), testu.begin() + u.size());
//
//                form_dense_jacobian_fd(fespace, disc, geo_map, std::span<T>(testu), std::span<T>(testres), densejac);
//
//                for(IDX i = 0; i < res.size() + mdg_res.size(); ++i){
//                    for(IDX j = 0; j < u.size() + coord.size(); ++j){
//                        T petsc_val;
//                        MatGetValue(jac, i, j, &petsc_val);
//                        if(std::abs( densejac[i, j] - petsc_val ) > 1e-6){
//                            std::cout << "i: " << i << " | j: " << j << " mismatch -- petsc: " << petsc_val 
//                                << " | dense:: " << densejac[i, j] << std::endl;
//                        }
//                    }
//                }
//                
            }

            // assume nonzero structure of jacobian remains unchanged
            if(explicitly_form_subproblem){
                MatProductSymbolic(subproblem_mat);
            }

            // set the initial residual norm
            PetscCallAbort(PETSC_COMM_WORLD, VecNorm(res_data, NORM_2, &(conv_criteria.r0)));

            IDX k;
            for(k = 0; k < conv_criteria.kmax; ++k){


//                PetscViewer jacobian_viewer;
//                PetscViewerASCIIOpen(PETSC_COMM_WORLD, ("iceicle_data/jacobian_view" + std::to_string(k) + ".dat").c_str(), &jacobian_viewer);
//                PetscViewerPushFormat(jacobian_viewer, PETSC_VIEWER_ASCII_DENSE);
//                MatView(jac, jacobian_viewer);
//
//                PetscViewer r_viewer;
//                PetscViewerASCIIOpen(PETSC_COMM_WORLD, ("iceicle_data/residual_view" + std::to_string(k) + ".dat").c_str(), &r_viewer);
//                PetscViewerPushFormat(r_viewer, PETSC_VIEWER_ASCII_DENSE);
//                VecView(res_data, r_viewer);


                // Form the subproblem
                if(explicitly_form_subproblem){

                    // JTJ 
                    MatProductNumeric(subproblem_mat);

                    // Regularization
                    Vec lambda;
                    VecCreate(PETSC_COMM_WORLD, &lambda);
                    VecSetSizes(lambda, u_layout.size() + geo_layout.size(), PETSC_DETERMINE);
                    VecSetFromOptions(lambda);
                    {
                        petsc::VecSpan lambda_view{lambda};
                        PetscCallAbort(PETSC_COMM_WORLD,
                                MatGetColumnNorms(jac, NORM_2, lambda_view.data()));

                        for(PetscInt i = 0; i < u_layout.size(); ++i)
                            { lambda_view[i] *= lambda_u; }
                        for(PetscInt i = u_layout.size(); i < u_layout.size() + geo_layout.size(); ++i)
                            { lambda_view[i] = std::max(lambda_b, lambda_view[i] * lambda_b); }

                        for(auto el : fespace.elements){

                            // get the minimum jacobian determinant of all quadrature points
                            T detJ = 1.0;
                            for(int igauss = 0; igauss < el.nQP(); ++igauss){
                                auto J = el.geo_el->Jacobian(fespace.meshptr->nodes, el.getQP(igauss).abscisse);
                                detJ = std::min(std::abs(detJ), std::abs(NUMTOOL::TENSOR::FIXED_SIZE::determinant(J)));

                            }
                            detJ = std::max(1e-8, std::abs(detJ));
                            IDX nodes_size = geo_layout.geo_map.selected_nodes.size();
                            for(auto inode : el.geo_el->nodes_span()){
                                IDX geo_dof = geo_layout.geo_map.inv_selected_nodes[inode];
                                if(geo_dof != nodes_size){
                                    for(int iv = 0; iv < geo_layout.nv(geo_dof); ++iv){
                                        lambda_view[u_layout.size() + geo_layout[geo_dof, iv]] += lambda_lag / (detJ);
                                    }
                                }
                            }
                        }
                    }
                    MatDiagonalSet(subproblem_mat, lambda, ADD_VALUES);
                    VecDestroy(&lambda);
                } else {
                    // TODO: 
                }
                PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyBegin(subproblem_mat, MAT_FINAL_ASSEMBLY));
                PetscCallAbort(PETSC_COMM_WORLD, MatAssemblyEnd(subproblem_mat, MAT_FINAL_ASSEMBLY));

                // form JTr 
                PetscCallAbort(PETSC_COMM_WORLD, MatMultTranspose(jac, res_data, Jtr));

                // Solve the subproblem
                PetscCallAbort(PETSC_COMM_WORLD, KSPSetOperators(ksp, subproblem_mat, subproblem_mat));
                PetscCallAbort(PETSC_COMM_WORLD, KSPSolve(ksp, Jtr, du_data));

                // update u 
                if constexpr (std::is_same_v<ls_type, no_linesearch<T, IDX>>){
                    petsc::VecSpan du_view{du_data};
                    fespan du{du_view.data(), u.get_layout()};
                    axpy(-1.0, du, u);

                    // x update
                    component_span dx{du_view.data() + u.size(), geo_layout};
                    axpy(-1.0, dx, coord);
                } else {
                    // its linesearchin time!
                    
                    T alpha = linesearch([&](T alpha_arg){

                        // === Set up working arrays for linesearch update ===

                        // u step for linesearch
                        std::vector<T> u_step_storage(u.size());
                        fespan u_step{u_step_storage.data(), u.get_layout()};
                        copy_fespan(u, u_step);

                        // x step for linesearch
                        std::vector<T> x_step_storage(geo_layout.size());
                        component_span x_step{x_step_storage, geo_layout};

                        // working array for linesearch residuals
                        std::vector<T> r_work_storage(u.size());
                        fespan res_work{r_work_storage.data(), u.get_layout()};

                        std::vector<T> r_mdg_work_storage(ic_layout.size());
                        dofspan mdg_res{r_mdg_work_storage, ic_layout};

                        // === Compute the Update ===

                        petsc::VecSpan du_view{du_data};
                        fespan du{du_view.data(), u.get_layout()};
                        axpy(-alpha_arg, du, u_step);

                        // x update
                        component_span dx{du_view.data() + u.size(), geo_layout};
                        axpy(-alpha_arg, dx, coord);

                        // === Get the residuals ===
                        form_residual(fespace, disc, u_step, res_work);
                        form_mdg_residual(fespace, disc, u_step, geo_map, mdg_res);
                        T rnorm = res_work.vector_norm() + mdg_res.vector_norm();

                        // Add any penalties

                        return rnorm;
                    });


                    // Perform the linesearch update
                    petsc::VecSpan du_view{du_data};
                    fespan du{du_view.data(), u.get_layout()};
                    axpy(-alpha, du, u);

                    // x update
                    component_span dx{du_view.data() + u.size(), geo_layout};
                    axpy(-alpha, dx, coord);
                }

                // clear out matrices 
                MatZeroEntries(jac); // zero out the jacobian

                // get updated residual and jacobian 
                {
                    petsc::VecSpan res_view{res_data};
                    fespan res{res_view.data(), u.get_layout()};
                    form_petsc_jacobian_fd(fespace, disc, u, res, jac);
                    dofspan mdg_res{res_view.data() + u_layout.size(), ic_layout};
                    form_petsc_mdg_jacobian_fd(fespace, disc, u, coord, mdg_res, jac);

                    // assemble the Jacobian matrix 
                    MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
                    MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);

                    // Test Jacobian correctness
//                    std::vector<T> densejac_data( (res.size() + mdg_res.size()) * (u.size() + coord.size()) );
//                    std::mdspan<T, std::extents<std::size_t, std::dynamic_extent, std::dynamic_extent>, std::layout_right> densejac{
//                        densejac_data.data(), res.size() + mdg_res.size(), u.size() + coord.size()
//                    };
//                    std::vector<T> testres(res.size() + mdg_res.size());
//                    std::vector<T> testu(u.size() + coord.size());
//                    std::copy_n(u.data(), u.size(), testu.begin());
//                    std::copy_n(coord.data(), coord.size(), testu.begin() + u.size());
//
//                    form_dense_jacobian_fd(fespace, disc, geo_map, std::span<T>(testu), std::span<T>(testres), densejac);
//
//                    for(IDX i = 0; i < res.size() + mdg_res.size(); ++i){
//                        for(IDX j = 0; j < u.size() + coord.size(); ++j){
//                            T petsc_val;
//                            MatGetValue(jac, i, j, &petsc_val);
//                            if(std::abs( densejac[i, j] - petsc_val ) > 1e-6 ){ 
//                                std::cout << "i: " << i << " | j: " << j << " mismatch -- petsc: " << petsc_val 
//                                    << " | dense:: " << densejac[i, j] << std::endl;
//                            }
//                        }
//                    }
                }

                // get the residual norm
                T rk;
                PetscCallAbort(PETSC_COMM_WORLD, VecNorm(res_data, NORM_2, &rk));

                // Diagnostics 
                if(idiag > 0 && k % idiag == 0) {
                    diag_callback(k, res_data, du_data);
                }

                // visualization
                if(ivis > 0 && k % ivis == 0) {
                    vis_callback(k, res_data, du_data);
                }

                // test convergence
                if(conv_criteria.done_callback(rk)) break;
            }
            return k;
        }
    };

}
