/**
 * @brief A finite element space represents a collection of finite elements 
 * and trace spaces (if applicable)
 * that provides a general interface to a finite element discretization of the domain
 * and simple generation utilities
 *
 * @author Gianni Absillis (gabsill@ncsu.edu)
 */
#pragma once

#include "iceicle/basis/basis.hpp"
#include "iceicle/crs.hpp"
#include "iceicle/element/finite_element.hpp"
#include <iceicle/element/reference_element.hpp>
#include "iceicle/fe_definitions.hpp"
#include "iceicle/fe_function/dglayout.hpp"
#include "iceicle/fe_function/cg_map.hpp"
#include "iceicle/geometry/face.hpp"
#include "iceicle/geometry/geo_element.hpp"
#include "iceicle/quadrature/QuadratureRule.hpp"
#include <iceicle/mesh/mesh.hpp>
#include <iceicle/tmp_utils.hpp>
#include <Numtool/tmp_flow_control.hpp>
#include <map>
#include <type_traits>

#ifdef ICEICLE_USE_MPI 
#include <mpi.h>
#endif

namespace iceicle {
    /**
     * Key to define the surjective mapping from an element 
     * to the corresponding evaluation
     */
    struct FETypeKey {
        DOMAIN_TYPE domain_type;

        int basis_order;

        int geometry_order;

        FESPACE_ENUMS::FESPACE_QUADRATURE qtype;

        FESPACE_ENUMS::FESPACE_BASIS_TYPE btype;


        friend bool operator<(const FETypeKey &l, const FETypeKey &r){
            using namespace FESPACE_ENUMS;
            if(l.qtype != r.qtype){
                return (int) l.qtype < (int) r.qtype;
            } else if(l.btype != r.btype){
                return (int) l.btype < (int) r.btype;
            } else if(l.domain_type != r.domain_type) {
                return (int) l.domain_type < (int) r.domain_type;
            } else if (l.geometry_order != r.geometry_order){
                return l.geometry_order < r.geometry_order;
            } else if( l.basis_order != r.basis_order) {
                return l.basis_order < r.basis_order;
            } else {
                // they are equal so less than is false for both cases
                return false;
            }
        }
    };

    /**
     * key to define surjective mapping from trace space to 
     * corresponding evaluation 
     */
    struct TraceTypeKey {

        FESPACE_ENUMS::FESPACE_BASIS_TYPE btype_l;

        FESPACE_ENUMS::FESPACE_BASIS_TYPE btype_r;

        int basis_order_l;

        int basis_order_r;

        int basis_order_trace;

        int geometry_order;

        DOMAIN_TYPE domain_type;

        FESPACE_ENUMS::FESPACE_QUADRATURE qtype;

        unsigned int face_info_l;

        unsigned int face_info_r;

        auto operator <=>(const TraceTypeKey&) const = default;
    };

    /**
     * @brief Collection of FiniteElements and TraceSpaces
     * to form a unified interface for a discretization of a domain 
     *
     * @tparam T the numeric type 
     * @tparam IDX the index type 
     * @tparam ndim the number of dimensions
     */
    template<typename T, typename IDX, int ndim>
    class FESpace {
        public:

        using ElementType = FiniteElement<T, IDX, ndim>;
        using TraceType = TraceSpace<T, IDX, ndim>;
        using GeoElementType = GeometricElement<T, IDX, ndim>;
        using GeoFaceType = Face<T, IDX, ndim>;
        using MeshType = AbstractMesh<T, IDX, ndim>;
        using BasisType = Basis<T, ndim>;
        using QuadratureType = QuadratureRule<T, IDX, ndim>;

        /// @brief what type of finite element space is being represented
        enum class SPACE_TYPE {
            L2, /// L2 elements are fully discontinuous at the interfaces
            ISOPARAMETRIC_H1, /// A continuous finite element space over the whole domain 
                              /// Basis functions for solution are equivalent
                              /// to basis functionss for the geometry
        };

        /// @brief what type of finite element space is being represented
        /// by this FESpace
        SPACE_TYPE type;

        /// @brief pointer to the mesh used
        MeshType *meshptr;

        /// @brief Array of finite elements in the space
        std::vector<ElementType> elements;

        /// @brief Array of trace spaces in the space 
        std::vector<TraceType> traces;

        /// @brief the start index of the interior traces 
        std::size_t interior_trace_start;
        /// @brief the end index of the interior traces (exclusive) 
        std::size_t interior_trace_end;

        /// @brief the start index of the boundary traces 
        std::size_t bdy_trace_start;
        /// @brief the end index of the boundary traces (exclusive)
        std::size_t bdy_trace_end;

        /** @brief maps local dofs to global dofs for dg space */
        dg_dof_map<IDX> dg_map;

        /** @brief maps local dofs to global dofs for cg space */
        cg_dof_map<T, IDX, ndim> cg_map;

        /** @brief the mapping of faces connected to each node */
        util::crs<IDX> fac_surr_nodes;

        /** @brief the mapping of elements connected to each node */
        util::crs<IDX, IDX> el_surr_nodes;

        /** @brief the mapping of faces connected to each element */
        util::crs<IDX> fac_surr_el;

        /// @brief element information recieved from each respective MPI rank
        std::vector<std::vector<ElementType>> comm_elements;

        private:

        // ========================================
        // = Maps to Basis, Quadrature, and Evals =
        // ========================================

        using ReferenceElementType = ReferenceElement<T, IDX, ndim>;
        using ReferenceTraceType = ReferenceTraceSpace<T, IDX, ndim>;
        std::map<FETypeKey, ReferenceElementType> ref_el_map;
        std::map<TraceTypeKey, ReferenceTraceType> ref_trace_map;

        public:

        // default constructor
        FESpace() = default;

        // delete copy semantics
        FESpace(const FESpace &other) = delete;
        FESpace<T, IDX, ndim>& operator=(const FESpace &other) = delete;

        // keep move semantics
        FESpace(FESpace &&other) = default;
        FESpace<T, IDX, ndim>& operator=(FESpace &&other) = default;

        /**
         * @brief construct an FESpace with uniform 
         * quadrature rules, and basis functions over all elements 
         *
         * @tparam basis_order the polynomial order of 1D basis functions
         *
         * @param meshptr pointer to the mesh 
         * @param basis_type enumeration of what basis to use 
         * @param quadrature_type enumeration of what quadrature rule to use 
         * @param basis_order_arg for template argument deduction of the basis order
         */
        template<int basis_order>
        FESpace(
            MeshType *meshptr,
            FESPACE_ENUMS::FESPACE_BASIS_TYPE basis_type,
            FESPACE_ENUMS::FESPACE_QUADRATURE quadrature_type,
            tmp::compile_int<basis_order> basis_order_arg
        ) : type{SPACE_TYPE::L2}, meshptr(meshptr), cg_map{*meshptr}, elements{} {

            // Generate the Finite Elements
            elements.reserve(meshptr->nelem());
            for(ElementTransformation<T, IDX, ndim>* geo_trans : meshptr->el_transformations){
                // create the Element Domain type key
                FETypeKey fe_key = {
                    .domain_type = geo_trans->domain_type,
                    .basis_order = basis_order,
                    .geometry_order = geo_trans->order,
                    .qtype = quadrature_type,
                    .btype = basis_type
                };

                // check if an evaluation doesn't exist yet
                if(ref_el_map.find(fe_key) == ref_el_map.end()){
                    ref_el_map[fe_key] = ReferenceElementType(geo_trans->domain_type, geo_trans->order, basis_type, quadrature_type, basis_order_arg);
                }
                ReferenceElementType &ref_el = ref_el_map[fe_key];

                // this will be the index of the new element
                IDX ielem = elements.size();

                // create the finite element
                ElementType fe{
                    .trans = geo_trans, 
                    .basis = ref_el.basis.get(),
                    .quadrule = ref_el.quadrule.get(),
                    .qp_evals = std::span<const BasisEvaluation<T, ndim>>{ref_el.evals},
                    .inodes = meshptr->conn_el.rowspan(ielem), // NOTE: meshptr cannot invalidate anymore
                    .coord_el = meshptr->coord_els.rowspan(ielem),
                    .elidx = ielem
                };

                // add to the elements list
                elements.push_back(fe);
            }


#ifdef ICEICLE_USE_MPI
            // ========================
            // = Communicate Elements =
            // ========================
            int myrank, nrank;
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
            MPI_Comm_size(MPI_COMM_WORLD, &nrank);

            comm_elements.resize(nrank);

            // set up the communicated FiniteElements only if more than 1 process
            for(int irank = 0; irank < nrank; ++irank){
                // basis order, quadrature_type, and basis_type we know
                int geometry_order, domain_type;
                for(auto& geo_el_info : meshptr->communicated_elements[irank]){
                    
                    // create the Element Domain type key
                    FETypeKey fe_key = {
                        .domain_type = geo_el_info.trans->domain_type,
                        .basis_order = basis_order,
                        .geometry_order = geo_el_info.trans->order,
                        .qtype = quadrature_type,
                        .btype = basis_type
                    };

                    // check if an evaluation doesn't exist yet
                    if(ref_el_map.find(fe_key) == ref_el_map.end()){
                        ref_el_map[fe_key] = ReferenceElementType(geo_el_info.trans->domain_type,
                                geo_el_info.trans->order, basis_type, quadrature_type, basis_order_arg);
                    }
                    ReferenceElementType &ref_el = ref_el_map[fe_key];
                
                    // create the finite element
                    ElementType fe(
                        geo_el_info.trans,
                        ref_el.basis.get(),
                        ref_el.quadrule.get(),
                        std::span<const BasisEvaluation<T, ndim>>{ref_el.evals},
                        geo_el_info.conn_el,
                        geo_el_info.coord_el,
                        elements.size() // this will be the index of the new element
                    );

                    comm_elements[irank].push_back(fe);
                }
            }
#endif

            // Generate the Trace Spaces
            traces.reserve(meshptr->faces.size());
            for(const auto& fac : meshptr->faces){
                // NOTE: assuming element indexing is the same as the mesh still

                bool is_interior = fac->bctype == BOUNDARY_CONDITIONS::INTERIOR;
                ElementType *elptrL = &elements[fac->elemL];
                ElementType *elptrR = (is_interior) ? &elements[fac->elemR] : &elements[fac->elemL];

#ifdef ICEICLE_USE_MPI
            // update the boundary trace spaces to use the new communicated FiniteElements
                if(fac->bctype == BOUNDARY_CONDITIONS::PARALLEL_COM) {
                    auto [jrank, imleft] = decode_mpi_bcflag(fac->bcflag);
                    IDX jlocal_elidx = (imleft) ? fac->elemR : fac->elemL;

                    std::vector<IDX> &comm_el_idxs = meshptr->el_recv_list[jrank];
                    auto itr = lower_bound(comm_el_idxs.begin(), comm_el_idxs.end(), jlocal_elidx);
                    std::size_t index = distance(comm_el_idxs.begin(), itr);

                    if(imleft){
                        elptrR = &(comm_elements[jrank].at(index));
                    } else {
                        elptrL = &(comm_elements[jrank].at(index));
                        elptrR = &elements[fac->elemR]; // special case because parallel faces are essential interior
                    }
                }
#endif
                ElementType& elL = *elptrL;
                ElementType& elR = *elptrR;

                int geo_order = std::max(elL.trans->order, elR.trans->order);

                auto geo_order_dispatch = [&]<int geo_order>() -> int{
                    TraceTypeKey trace_key = { 
                        .basis_order_l = elL.basis->getPolynomialOrder(),
                        .basis_order_r = elR.basis->getPolynomialOrder(),
                        .basis_order_trace = std::max(elL.basis->getPolynomialOrder(), elR.basis->getPolynomialOrder()), 
                        .geometry_order = geo_order,
                        .domain_type = fac->domain_type(),
                        .qtype = quadrature_type,
                        .face_info_l = fac->face_infoL,
                        .face_info_r = fac->face_infoR
                    };

                    if(ref_trace_map.find(trace_key) == ref_trace_map.end()){
                        ref_trace_map[trace_key] = ReferenceTraceType(fac.get(),
                            basis_type, quadrature_type, 
                            *(elL.basis), *(elR.basis),
                            std::integral_constant<int, basis_order>{},
                            std::integral_constant<int, geo_order>{});
                    }
                    ReferenceTraceType &ref_trace = ref_trace_map[trace_key];
                    
                    if(is_interior || fac->bctype == BOUNDARY_CONDITIONS::PARALLEL_COM){
                        // parallel bdy faces are essentially also interior faces 
                        // aside from being a bit *special* :3
                        TraceType trace{ fac.get(), &elL, &elR, ref_trace.trace_basis.get(),
                            ref_trace.quadrule.get(),
                            std::span<const BasisEvaluation<T, ndim>>{ref_trace.evals_l},
                            std::span<const BasisEvaluation<T, ndim>>{ref_trace.evals_r},
                            (IDX) traces.size() };
                        traces.push_back(trace);
                    } else {
                        TraceType trace = TraceType::make_bdy_trace_space(
                            fac.get(), &elL, ref_trace.trace_basis.get(), 
                            ref_trace.quadrule.get(), 
                            std::span<const BasisEvaluation<T, ndim>>{ref_trace.evals_l},
                            std::span<const BasisEvaluation<T, ndim>>{ref_trace.evals_r},
                            (IDX) traces.size());
                        traces.push_back(trace);
                    }

                    return 0;
                };

                NUMTOOL::TMP::invoke_at_index(
                    NUMTOOL::TMP::make_range_sequence<int, 1, MAX_DYNAMIC_ORDER>{},
                    geo_order,
                    geo_order_dispatch                    
                );
            }
            // reuse the face indexing from the mesh
            interior_trace_start = meshptr->interiorFaceStart;
            interior_trace_end = meshptr->interiorFaceEnd;
            bdy_trace_start = meshptr->bdyFaceStart;
            bdy_trace_end = meshptr->bdyFaceEnd;

            // generate the dof offsets 
            dg_map = dg_dof_map{elements};

            // ===================================
            // = Build the connectivity matrices =
            // ===================================

            // generate the face surrounding nodes connectivity matrix 
            std::vector<std::vector<IDX>> connectivity_ragged(meshptr->n_nodes());
            for(int itrace = 0; itrace < traces.size(); ++itrace){
                const TraceType& trace = traces[itrace];
                for(IDX inode : trace.face->nodes_span()){
                    connectivity_ragged[inode].push_back(itrace);
                }
            }
            fac_surr_nodes = util::crs{connectivity_ragged};

            el_surr_nodes = util::crs{meshptr->elsup};

            std::vector<std::vector<IDX>> fac_surr_el_ragged(elements.size());
            for(int itrace = 0; itrace < traces.size(); ++itrace) {
                const TraceType& trace = traces[itrace];
                if(trace.face->bctype == BOUNDARY_CONDITIONS::PARALLEL_COM){
                    // take some extra care to not add the wrong element index
                    auto [jrank, imleft] = decode_mpi_bcflag(trace.face->bcflag);
                    if(imleft){
                        fac_surr_el_ragged[trace.elL.elidx].push_back(itrace);
                    } else {
                        fac_surr_el_ragged[trace.elR.elidx].push_back(itrace);
                    }
                } else {
                    fac_surr_el_ragged[trace.elL.elidx].push_back(itrace);
                    fac_surr_el_ragged[trace.elR.elidx].push_back(itrace);
                }
            }
            fac_surr_el = util::crs{fac_surr_el_ragged};
        } 

        /// @brief construct an FESpace that represents an isoparametric CG space
        /// to the given mesh 
        /// @param meshptr pointer to the mesh
        FESpace(MeshType *meshptr) : type{SPACE_TYPE::ISOPARAMETRIC_H1}, meshptr(meshptr), cg_map{*meshptr}, elements{} {
            
            // Generate the Finite Elements
            elements.reserve(meshptr->nelem());
            for(ElementTransformation<T, IDX, ndim>* geo_trans : meshptr->el_transformations){
                // create the Element Domain type key
                FETypeKey fe_key = {
                    .domain_type = geo_trans->domain_type,
                    .basis_order = geo_trans->order,
                    .geometry_order = geo_trans->order,
                    .qtype = FESPACE_ENUMS::FESPACE_QUADRATURE::GAUSS_LEGENDRE,
                    .btype = FESPACE_ENUMS::FESPACE_BASIS_TYPE::LAGRANGE 
                };

                // check if an evaluation doesn't exist yet
                if(ref_el_map.find(fe_key) == ref_el_map.end()){
                    ref_el_map[fe_key] = ReferenceElementType(geo_trans->domain_type, geo_trans->order);
                }
                ReferenceElementType &ref_el = ref_el_map[fe_key];
               
                // this will be the index of the new element
                IDX ielem = elements.size();

                // create the finite element
                ElementType fe{
                    .trans = geo_trans, 
                    .basis = ref_el.basis.get(),
                    .quadrule = ref_el.quadrule.get(),
                    .qp_evals = std::span<const BasisEvaluation<T, ndim>>{ref_el.evals},
                    .inodes = meshptr->conn_el.rowspan(ielem), // NOTE: meshptr cannot invalidate anymore
                    .coord_el = meshptr->coord_els.rowspan(ielem),
                    .elidx = ielem
                };

                // add to the elements list
                elements.push_back(fe);
            }

#ifdef ICEICLE_USE_MPI
            // ========================
            // = Communicate Elements =
            // ========================
            int myrank, nrank;
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
            MPI_Comm_size(MPI_COMM_WORLD, &nrank);

            comm_elements.resize(nrank);

            // set up the communicated FiniteElements only if more than 1 process
            for(int irank = 0; irank < nrank; ++irank){
                // basis order, quadrature_type, and basis_type we know
                int geometry_order, domain_type;
                for(auto& comm_el : meshptr->communicated_elements[irank]){
                    
                    FETypeKey fe_key = {
                        .domain_type = comm_el.trans->domain_type,
                        .basis_order = comm_el.trans->order,
                        .geometry_order = comm_el.trans->order,
                        .qtype = FESPACE_ENUMS::FESPACE_QUADRATURE::GAUSS_LEGENDRE,
                        .btype = FESPACE_ENUMS::FESPACE_BASIS_TYPE::LAGRANGE 
                    };

                    // check if an evaluation doesn't exist yet
                    if(ref_el_map.find(fe_key) == ref_el_map.end()){
                        ref_el_map[fe_key] = ReferenceElementType(comm_el.trans->domain_type, comm_el.trans->order);
                    }
                    ReferenceElementType &ref_el = ref_el_map[fe_key];
                
                    // this will be the index of the new element
                    IDX ielem = elements.size();

                    // create the finite element
                    ElementType fe{
                        .trans = comm_el.trans, 
                        .basis = ref_el.basis.get(),
                        .quadrule = ref_el.quadrule.get(),
                        .qp_evals = std::span<const BasisEvaluation<T, ndim>>{ref_el.evals},
                        .inodes = comm_el.conn_el, // NOTE: meshptr cannot invalidate anymore
                        .coord_el = comm_el.coord_el,
                        .elidx = ielem
                    };

                    comm_elements[irank].push_back(fe);
                }
            }
#endif
            // Generate the Trace Spaces
            traces.reserve(meshptr->faces.size());
            for(const auto& fac : meshptr->faces){
                // NOTE: assuming element indexing is the same as the mesh still

                bool is_interior = fac->bctype == BOUNDARY_CONDITIONS::INTERIOR;
                ElementType *elptrL = &elements[fac->elemL];
                ElementType *elptrR = (is_interior) ? &elements[fac->elemR] : &elements[fac->elemL];

#ifdef ICEICLE_USE_MPI
            // update the boundary trace spaces to use the new communicated FiniteElements
                if(fac->bctype == BOUNDARY_CONDITIONS::PARALLEL_COM) {
                    auto [jrank, imleft] = decode_mpi_bcflag(fac->bcflag);
                    IDX jlocal_elidx = (imleft) ? fac->elemR : fac->elemL;

                    std::vector<IDX> &comm_el_idxs = meshptr->el_recv_list[jrank];
                    auto itr = lower_bound(comm_el_idxs.begin(), comm_el_idxs.end(), jlocal_elidx);
                    std::size_t index = distance(comm_el_idxs.begin(), itr);

                    if(imleft){
                        elptrR = &(comm_elements[jrank].at(index));
                    } else {
                        elptrL = &(comm_elements[jrank].at(index));
                        elptrR = &elements[fac->elemR]; // special case because parallel faces are essential interior
                    }
                }
#endif
                ElementType& elL = *elptrL;
                ElementType& elR = *elptrR;

                int geo_order = std::max(elL.trans->order, elR.trans->order);

                auto geo_order_dispatch = [&]<int geo_order>() -> int{
                    TraceTypeKey trace_key = { 
                        .basis_order_l = elL.basis->getPolynomialOrder(),
                        .basis_order_r = elR.basis->getPolynomialOrder(),
                        .basis_order_trace = std::max(elL.basis->getPolynomialOrder(), elR.basis->getPolynomialOrder()), 
                        .geometry_order = geo_order,
                        .domain_type = fac->domain_type(),
                        .qtype = FESPACE_ENUMS::FESPACE_QUADRATURE::GAUSS_LEGENDRE,
                        .face_info_l = fac->face_infoL,
                        .face_info_r = fac->face_infoR
                    };

                    if(ref_trace_map.find(trace_key) == ref_trace_map.end()){
                        ref_trace_map[trace_key] = ReferenceTraceType(
                            fac.get(),
                            FESPACE_ENUMS::FESPACE_BASIS_TYPE::LAGRANGE,
                            FESPACE_ENUMS::FESPACE_QUADRATURE::GAUSS_LEGENDRE,
                            *(elL.basis),
                            *(elR.basis),
                            std::integral_constant<int, geo_order>{},
                            std::integral_constant<int, geo_order>{});
                    }
                    ReferenceTraceType &ref_trace = ref_trace_map[trace_key];
                    
                    if(is_interior || fac->bctype == BOUNDARY_CONDITIONS::PARALLEL_COM){
                        // parallel bdy faces are essentially also interior faces 
                        // aside from being a bit *special* :3
                        TraceType trace{ fac.get(), &elL, &elR, ref_trace.trace_basis.get(),
                            ref_trace.quadrule.get(), 
                            std::span<const BasisEvaluation<T, ndim>>{ref_trace.evals_l},
                            std::span<const BasisEvaluation<T, ndim>>{ref_trace.evals_r},
                            (IDX) traces.size() };
                        traces.push_back(trace);
                    } else {
                        TraceType trace = TraceType::make_bdy_trace_space(
                            fac.get(), &elL, ref_trace.trace_basis.get(), 
                            ref_trace.quadrule.get(),
                            std::span<const BasisEvaluation<T, ndim>>{ref_trace.evals_l},
                            std::span<const BasisEvaluation<T, ndim>>{ref_trace.evals_r},
                            (IDX) traces.size());
                        traces.push_back(trace);
                    }

                    return 0;
                };

                NUMTOOL::TMP::invoke_at_index(
                    NUMTOOL::TMP::make_range_sequence<int, 1, MAX_DYNAMIC_ORDER>{},
                    geo_order,
                    geo_order_dispatch                    
                );
            }
            // reuse the face indexing from the mesh
            interior_trace_start = meshptr->interiorFaceStart;
            interior_trace_end = meshptr->interiorFaceEnd;
            bdy_trace_start = meshptr->bdyFaceStart;
            bdy_trace_end = meshptr->bdyFaceEnd;

            // generate the dof offsets 
            dg_map = dg_dof_map{elements};

            // ===================================
            // = Build the connectivity matrices =
            // ===================================

            // generate the face surrounding nodes connectivity matrix 
            std::vector<std::vector<IDX>> connectivity_ragged(meshptr->n_nodes());
            for(int itrace = 0; itrace < traces.size(); ++itrace){
                const TraceType& trace = traces[itrace];
                for(IDX inode : trace.face->nodes_span()){
                    connectivity_ragged[inode].push_back(itrace);
                }
            }
            fac_surr_nodes = util::crs{connectivity_ragged};

            el_surr_nodes = util::crs{meshptr->elsup};

            std::vector<std::vector<IDX>> fac_surr_el_ragged(elements.size());
            for(int itrace = 0; itrace < traces.size(); ++itrace) {
                const TraceType& trace = traces[itrace];
                if(trace.face->bctype == BOUNDARY_CONDITIONS::PARALLEL_COM){
                    // take some extra care to not add the wrong element index
                    auto [jrank, imleft] = decode_mpi_bcflag(trace.face->bcflag);
                    if(imleft){
                        fac_surr_el_ragged[trace.elL.elidx].push_back(itrace);
                    } else {
                        fac_surr_el_ragged[trace.elR.elidx].push_back(itrace);
                    }
                } else {
                    fac_surr_el_ragged[trace.elL.elidx].push_back(itrace);
                    fac_surr_el_ragged[trace.elR.elidx].push_back(itrace);
                }
            }
            fac_surr_el = util::crs{fac_surr_el_ragged};
        }

        /**
         * @brief get the number of dg degrees of freedom in the entire fespace 
         * multiply this by the nummber of components to get the size requirement for 
         * a dg fespan or use the built_in function in the dg_map member
         * @return the number of dg degrees of freedom
         */
        constexpr std::size_t ndof_dg() const noexcept
        {
            return dg_map.calculate_size_requirement(1);
        }

        /**
         * @brief get the span that is the subset of the trace space list 
         * that only includes interior traces 
         * @return span over the interior traces 
         */
        std::span<TraceType> get_interior_traces(){
            return std::span<TraceType>{traces.begin() + interior_trace_start,
                traces.begin() + interior_trace_end};
        }

        /**
         * @brief get the span that is the subset of the trace space list 
         * that only includes boundary traces 
         * @return span over the boundary traces 
         */
        std::span<TraceType> get_boundary_traces(){
            return std::span<TraceType>{traces.begin() + bdy_trace_start,
                traces.begin() + bdy_trace_end};
        }

        auto print_info(std::ostream& out)
        -> std::ostream& {
            out << "Finite Element Space" << std::endl;
            switch(type){
                case SPACE_TYPE::L2:
                    out << "Space Type: ";
                    out << "L2" << std::endl;
                    out << "ndof: " << dg_map.size() << std::endl;
                    break;
                case SPACE_TYPE::ISOPARAMETRIC_H1:
                    out << "Space Type: ";
                    out << "H1 (isoparametric)" << std::endl;
                    out << "ndof: " << cg_map.size() << std::endl;
                    break;
            }
            return out;
        }

    };
    
}
