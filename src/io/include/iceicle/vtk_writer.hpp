#pragma once
#include "iceicle/anomaly_log.hpp"
#include "iceicle/fespace/fespace.hpp"
#include "iceicle/build_config.hpp"
#include "iceicle/fe_definitions.hpp"
#include "iceicle/iceicle_mpi_utils.hpp"
#include <iceicle/mesh/mesh.hpp>
#include <vtkLagrangeQuadrilateral.h>
#include <vtkLagrangeHexahedron.h>
#include <vtkLagrangeTriangle.h>
#include <vtkLagrangeCurve.h>
#include <vtkLagrangeTetra.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#include <filesystem>
namespace iceicle::io {

    /**
     * @brief a writer for pvtu files (vtk's parallel unstructured grid format)
     */
    template<typename T, typename IDX, int ndim, int conformity>
    class PVTUWriter {

        /// @brief abstract representation of a set of fields that can be written to a file
        /// and the mechanisms for writing it to a grid
        struct writeable_fieldset {
            
            struct FieldsetConcept{
                public:
                    virtual ~FieldsetConcept() = default;
                    virtual void do_write_data(const FESpace<T, IDX, ndim, conformity>& fespace, 
                            vtkSmartPointer<vtkUnstructuredGrid> vtk_grid) const = 0;
                    virtual auto clone() const -> std::unique_ptr<FieldsetConcept> = 0;
            };

            template< class FieldsetT >
            struct FieldsetModel : public FieldsetConcept {
                public:

                FieldsetT _fieldset;

                FieldsetModel(FieldsetT fieldset) : _fieldset{std::move(fieldset)}{}

                void do_write_data(const FESpace<T, IDX, ndim, conformity>& fespace, 
                        vtkSmartPointer<vtkUnstructuredGrid> vtk_grid) const override 
                { _fieldset.write_data(fespace, vtk_grid); }

                auto clone() const 
                -> std::unique_ptr<FieldsetConcept> override 
                { return std::make_unique<FieldsetModel>(*this); }

            };

            std::unique_ptr<FieldsetConcept> pimpl;

            public:

            writeable_fieldset() : pimpl{nullptr} {}

            template< typename FieldsetT >
            writeable_fieldset(FieldsetT fieldset)
            : pimpl{std::make_unique<FieldsetModel<FieldsetT>>(std::move(fieldset))}
            {}

            /// @brief check if the writer contains an actual writer value
            operator bool() const { return (bool) pimpl; }

            writeable_fieldset(const writeable_fieldset& other)
            : pimpl( other.pimpl->clone())
            {}

            writeable_fieldset& operator=(const writeable_fieldset& other) {
                // copy and swap
                writeable_fieldset tmp(other);
                std::swap(pimpl, tmp.pimpl);
                return *this;
            }

            // move 
            writeable_fieldset( writeable_fieldset&& other) = default;
            writeable_fieldset& operator=( writeable_fieldset&& other ) = default;

            void write_data(const FESpace<T, IDX, ndim, conformity>& fespace, 
                    vtkSmartPointer<vtkUnstructuredGrid> vtk_grid) const {
                if(pimpl)
                    pimpl->do_write_data(fespace, vtk_grid);
            }
        };

        struct mpi_rank_fieldset {
            void write_data(const FESpace<T, IDX, ndim, conformity>& fespace, 
                    vtkSmartPointer<vtkUnstructuredGrid> vtk_grid) const
            {
                // Create an Array to represent the MPI rank for each cell 
                vtkSmartPointer<vtkIntArray> rank_array = vtkSmartPointer<vtkIntArray>::New();
                rank_array->SetName("MPI rank");
                rank_array->SetNumberOfComponents(1);
                rank_array->SetNumberOfTuples(vtk_grid->GetCells()->GetNumberOfCells());
                for(IDX ielem = 0; ielem < fespace.elements.size(); ++ielem) {
                    rank_array->SetValue(ielem, mpi::mpi_world_rank());
                }
                vtk_grid->GetCellData()->AddArray(rank_array);
            }
        };

        // === Type aliases ===
        using Point = MATH::GEOMETRY::Point<T, ndim>;
        using FacePoint = MATH::GEOMETRY::Point<T, ndim - 1>;

        /// === VTK Reference Cells ===
        /// @brief the maximum polynomial order considered for output
        static constexpr int max_pn = std::max(build_config::FESPACE_BUILD_PN + 1,
                build_config::FESPACE_BUILD_GEO_PN + 1);
        std::vector< vtkSmartPointer< vtkCell > > reference_hypercubes;
        std::vector< vtkSmartPointer< vtkCell > > reference_simplices;
        std::vector< vtkSmartPointer< vtkCell > > reference_hypercube_faces;
        std::vector< vtkSmartPointer< vtkCell > > reference_simplex_faces;

        /// === VTK members ===

        /// @brief the vtk file writer 
        vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer;

        /// @brief the VTK mpi communicator
        vtkSmartPointer<vtkMPICommunicator> vtk_comm;

        // === General Members ===

        /// @brief the finite element space
        FESpace<T, IDX, ndim, conformity>& fespace;

        /// @brief all the fieldsets
        std::vector<writeable_fieldset> fieldsets;

        /// @brief the MPI communicator to use
        mpi::communicator_type comm;

        /// @brief the name of the data collection 
        std::string collection_name = "data";

        /// @brief the directory where the collection will be stored 
        std::filesystem::path data_directory;

        /// === Private implementation details ===

        /// @brief given an element transformation, the coordinates of the element, and solution order 
        /// generate the points to write to the vtk file for the given element
        [[nodiscard]] inline constexpr 
        auto get_ref_pts(
            const ElementTransformation<T, IDX, ndim>* trans, 
            const std::span<Point> el_coord,
            int order
        ) const
        -> std::pair<vtkSmartPointer<vtkCell>, std::vector<Point>>
        {
            int output_order = std::max(trans->order, order);
            switch(trans->domain_type){
                case DOMAIN_TYPE::HYPERCUBE:
                {
                    auto ref_cell = reference_hypercubes[output_order];
                    double* coords = ref_cell->GetParametricCoords();
                    std::vector<Point> refpts(ref_cell->GetNumberOfPoints());
                    for(int ipoin = 0; ipoin < ref_cell->GetNumberOfPoints(); ++ipoin){
                        double * coords_start = coords + 3 * ipoin;
                        MATH::GEOMETRY::Point<T, ndim> refpt{};
                        for(int idim = 0; idim < ndim; ++idim){
                            // shift from [0, 1] to [-1, 1]
                            refpt[idim] = coords_start[idim] * 2.0 - 1.0;
                        }
                        refpts[ipoin] = refpt;
                    }
                    return std::pair{ref_cell, refpts};
                } break;
                case DOMAIN_TYPE::SIMPLEX:
                {
                    auto ref_cell = reference_simplices[output_order];
                    double* coords = ref_cell->GetParametricCoords();
                    std::vector<Point> refpts(ref_cell->GetNumberOfPoints());
                    for(int ipoin = 0; ipoin < ref_cell->GetNumberOfPoints(); ++ipoin){
                        double * coords_start = coords + 3 * ipoin;
                        MATH::GEOMETRY::Point<T, ndim> refpt{};
                        refpts[ipoin] = refpt;
                    }
                    return std::pair{ref_cell, refpts};
                } break;
            }

            vtkSmartPointer<vtkCell> empty_cell = vtkSmartPointer<vtkLagrangeQuadrilateral>::New();
            util::AnomalyLog::log_anomaly("Could not get element points");
            return std::pair{empty_cell, std::vector<Point>{}};
        }

        /// @brief given a face and the solution order 
        /// generate the points to write to the vtk file for the given face 
        [[nodiscard]]inline
        auto get_ref_pts(const Face<T, IDX, ndim>* faceptr, int order) const
        -> std::pair<vtkSmartPointer<vtkCell>, std::vector<FacePoint>>
        {
            int output_order = std::max(faceptr->geometry_order(), order);
            switch(faceptr->domain_type()){
                case DOMAIN_TYPE::HYPERCUBE:
                {
                    auto ref_cell = reference_hypercube_faces[output_order];
                    double* coords = ref_cell->GetParametricCoords();
                    std::vector<FacePoint> refpts(ref_cell->GetNumberOfPoints());
                    for(int ipoin = 0; ipoin < ref_cell->GetNumberOfPoints(); ++ipoin){
                        double * coords_start = coords + 3 * ipoin;
                        FacePoint refpt{};
                        for(int idim = 0; idim < refpt.size(); ++idim){
                            // shift from [0, 1] to [-1, 1]
                            refpt[idim] = coords_start[idim] * 2.0 - 1.0;
                        }
                        refpts[ipoin] = refpt;
                    }
                    return std::pair{ref_cell, refpts};
                } break;
                case DOMAIN_TYPE::SIMPLEX:
                {
                    auto ref_cell = reference_simplex_faces[output_order];
                    double* coords = ref_cell->GetParametricCoords();
                    std::vector<FacePoint> refpts(ref_cell->GetNumberOfPoints());
                    for(int ipoin = 0; ipoin < ref_cell->GetNumberOfPoints(); ++ipoin){
                        double * coords_start = coords + 3 * ipoin;
                        FacePoint refpt{};
                        for(int idim = 0; idim < refpt.size(); ++idim){
                            refpt[idim] = coords_start[idim];
                        }
                        refpts[ipoin] = refpt;
                    }
                    return std::pair{ref_cell, refpts};
                } break;
            }

            util::AnomalyLog::log_anomaly("Could not get face points");
            vtkSmartPointer<vtkCell> empty_cell = vtkSmartPointer<vtkLagrangeQuadrilateral>::New();
            return std::pair{empty_cell, std::vector<FacePoint>{}};

        }

        /// @brief create and initialize a vtk grid 
        /// This will generate a grid with points for element-wise data, 
        /// and another section of points for trace-wise data
        ///
        /// Along with offsets to points for each element and trace, respectively
        [[nodiscard]] inline 
        auto create_vtk_grid() const
        -> std::tuple< vtkSmartPointer< vtkUnstructuredGrid >, std::vector<IDX>, std::vector<IDX> >
        {
            auto vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
            // allocate for the number of owned elements and all the traces
            vtk_grid->Allocate(fespace.elements.size() + fespace.traces.size()); 

            IDX ivtk_pt = 0;
            vtkNew<vtkPoints> vtk_coord;

            // ===================
            // = Write the Nodes =
            // ===================
            // This is done by looping through the elements and then traces,
            // getting the VTK reference cell,
            // then interpolating to vtk's reference domain points

            // generate points for each element
            std::vector<IDX> el_pts_offsets = {0};
            for(IDX ielem = 0; ielem < fespace.elements.size(); ++ielem){
                const FiniteElement<T, IDX, ndim>& el = fespace.elements[ielem];
                auto [ref_cell, ref_pts] = get_ref_pts(el.trans, 
                        el.coord_el, el.basis->getPolynomialOrder());

                std::vector<vtkIdType> pts{};
                for(Point ref_pt : ref_pts) {
                    Point act_pt = el.transform(ref_pt);
                    double vtk_point[3];
                    for(int idim = 0; idim < std::min(3, ndim); ++idim)
                        vtk_point[idim] = act_pt[idim];
                    for(int idim = ndim; idim < 3; ++idim)
                        vtk_point[idim] = 0.0;
                    vtk_coord->InsertPoint(ivtk_pt, vtk_point);
                    pts.push_back(ivtk_pt);
                    ++ivtk_pt;
                }
                el_pts_offsets.push_back(ivtk_pt);
                vtk_grid->InsertNextCell(ref_cell->GetCellType(),
                        ref_cell->GetNumberOfPoints(), pts.data());
            }

            // generate points for the faces
            std::vector<IDX> trace_pts_offsets = {ivtk_pt};
            for(IDX itrace = 0; itrace < fespace.traces.size(); ++itrace){
                const TraceSpace<T, IDX, ndim>& trace = fespace.traces[itrace];
                int order = std::max({trace.elL.basis->getPolynomialOrder(), trace.elR.basis->getPolynomialOrder(),
                        trace.trace_basis.getPolynomialOrder()});
                auto [ref_cell, ref_pts] = get_ref_pts(trace.face, order);

                std::vector<vtkIdType> pts{};
                for(FacePoint ref_pt : ref_pts) {
                    Point act_pt = trace.transform(ref_pt, fespace.meshptr->coord);
                    double vtk_point[3];
                    for(int idim = 0; idim < std::min(3, ndim - 1); ++idim)
                        vtk_point[idim] = act_pt[idim];
                    for(int idim = ndim; idim < 3; ++idim)
                        vtk_point[idim] = 0.0;
                    vtk_coord->InsertPoint(ivtk_pt, vtk_point);
                    pts.push_back(ivtk_pt);
                    ++ivtk_pt;
                }
                trace_pts_offsets.push_back(ivtk_pt);
                vtk_grid->InsertNextCell(ref_cell->GetCellType(),
                        ref_cell->GetNumberOfPoints(), pts.data());
            }
            vtk_grid->SetPoints(vtk_coord);

            return std::tuple{vtk_grid, el_pts_offsets, trace_pts_offsets};
        }


        public:

        // === Publically Accessible Type Aliases ===
        using value_type = T;

        // Constructor
        PVTUWriter(FESpace<T, IDX, ndim, conformity>& fespace, mpi::communicator_type comm)
        : writer{vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New()},
          reference_hypercubes(max_pn), reference_simplices(max_pn),
          reference_hypercube_faces(max_pn), reference_simplex_faces(max_pn),
          fespace{fespace}, fieldsets{mpi_rank_fieldset{}}, comm{comm}, 
          data_directory{std::filesystem::current_path() / "iceicle_data"}
        {
            // set up mpi stuffs
#ifdef ICEICLE_USE_MPI
            // set VTK to use our MPI communicator
            vtkSmartPointer<vtkMPICommunicator> vtk_comm = vtkSmartPointer<vtkMPICommunicator>::New();
            vtkMPICommunicatorOpaqueComm vtk_opaque_comm(&comm);
            vtk_comm->InitializeExternal(&vtk_opaque_comm);
            vtkSmartPointer<vtkMPIController> vtk_mpi_ctrl = vtkSmartPointer<vtkMPIController>::New();
            vtk_mpi_ctrl->SetCommunicator(vtk_comm);
            writer->SetController(vtk_mpi_ctrl);
#endif

            // set up reference elements
            for(int geo_order = 1; geo_order < max_pn; ++geo_order) {
                switch(ndim){
                    case 2:
                        {
                            auto quad = vtkSmartPointer<vtkLagrangeQuadrilateral>::New();
                            quad->SetOrder(geo_order, geo_order);
                            quad->Initialize();
                            reference_hypercubes[geo_order] = quad;
                            auto tri = vtkSmartPointer<vtkLagrangeTriangle>::New();
                            int npoin = (geo_order + 1) * (geo_order + 2) / 2;
                            tri->GetPointIds()->SetNumberOfIds(npoin);
                            tri->GetPoints()->SetNumberOfPoints(npoin);
                            tri->Initialize();
                            reference_simplices[geo_order] = tri;

                            auto segment = vtkSmartPointer<vtkLagrangeCurve>::New();
                            segment->GetPointIds()->SetNumberOfIds(geo_order + 1);
                            segment->GetPoints()->SetNumberOfPoints(geo_order + 1);
                            segment->Initialize();
                            reference_hypercube_faces[geo_order] = segment;
                            reference_simplex_faces[geo_order] = segment;
                        }
                        break;
                    case 3:
                        {
                            auto hex = vtkSmartPointer<vtkLagrangeHexahedron>::New();
                            hex->SetOrder(geo_order, geo_order, geo_order);
                            hex->Initialize();
                            reference_hypercubes[geo_order] = hex;
                            auto tetr = vtkSmartPointer<vtkLagrangeTetra>::New();
                            int npoin = (geo_order + 1) * (geo_order + 2) * (geo_order + 3) / 6;
                            tetr->GetPointIds()->SetNumberOfIds(npoin);
                            tetr->GetPoints()->SetNumberOfPoints(npoin);
                            tetr->Initialize();
                            reference_simplices[geo_order] = tetr;

                            { // face domains
                                auto quad = vtkSmartPointer<vtkLagrangeQuadrilateral>::New();
                                quad->SetOrder(geo_order, geo_order);
                                quad->Initialize();
                                reference_hypercube_faces[geo_order] = quad;
                                auto tri = vtkSmartPointer<vtkLagrangeTriangle>::New();
                                int npoin = (geo_order + 1) * (geo_order + 2) / 2;
                                tri->GetPointIds()->SetNumberOfIds(npoin);
                                tri->GetPoints()->SetNumberOfPoints(npoin);
                                tri->Initialize();
                                reference_simplex_faces[geo_order] = tri;
                            }
                        }
                        break;
                }
            }
        }

        /// @brief rename the collection
        void rename_collection(std::string_view new_name)
        { collection_name = new_name; }

        /// @brief write the mesh and field values in a .pvtu file
        /// @param itime the timestep 
        /// @param time the time value 
        /// NOTE: the user is responsible for making sure itim and time are unique
        void write(int itime, T time) const {

            auto [vtk_grid, element_pts_offsets, trace_pts_offsets] =
                create_vtk_grid();


            for(const writeable_fieldset& fieldset : fieldsets){
                fieldset.write_data(fespace, vtk_grid);
            }

            // write fields for time and cycle 
            vtkNew<vtkDoubleArray> t{};
            t->SetName("TIME");
            t->SetNumberOfTuples(1);
            t->SetTuple1(0, time);
            vtk_grid->GetFieldData()->AddArray(t);

            vtkNew<vtkDoubleArray> c{};
            c->SetName("CYCLE");
            c->SetNumberOfTuples(1);
            c->SetTuple1(0, itime);
            vtk_grid->GetFieldData()->AddArray(c);

            // write out to file
            writer->SetNumberOfPieces(mpi::mpi_world_size());
            writer->SetStartPiece(mpi::mpi_world_rank());
            writer->SetEndPiece(mpi::mpi_world_rank());

            writer->SetInputData(vtk_grid);
            writer->SetFileName(((data_directory / collection_name).string() 
                        + "_" + std::to_string(itime) + ".pvtu").c_str());
            writer->SetDataModeToAscii();
            writer->Write();
        }

    };

    template<class T, class IDX, int ndim>
    auto write_mesh_vtk(const AbstractMesh<T, IDX, ndim>& mesh) {

        int max_pn = build_config::FESPACE_BUILD_GEO_PN + 1;
        std::vector< vtkSmartPointer< vtkCell > >  reference_hypercubes(max_pn);
        std::vector< vtkSmartPointer< vtkCell > >  reference_simplices(max_pn);
        for(int geo_order = 1; geo_order < max_pn; ++geo_order) {
            switch(ndim){
                case 2:
                    {
                        auto quad = vtkSmartPointer<vtkLagrangeQuadrilateral>::New();
                        quad->SetOrder(geo_order, geo_order);
                        quad->Initialize();
                        reference_hypercubes[geo_order] = quad;
                        auto tri = vtkSmartPointer<vtkLagrangeTriangle>::New();
                        int npoin = (geo_order + 1) * (geo_order + 2) / 2;
                        tri->GetPointIds()->SetNumberOfIds(npoin);
                        tri->GetPoints()->SetNumberOfPoints(npoin);
                        tri->Initialize();
                        reference_simplices[geo_order] = tri;
                    }
                    break;
                case 3:
                    {
                        auto hex = vtkSmartPointer<vtkLagrangeHexahedron>::New();
                        hex->SetOrder(geo_order, geo_order, geo_order);
                        hex->Initialize();
                        reference_hypercubes[geo_order] = hex;
                        auto tetr = vtkSmartPointer<vtkLagrangeTetra>::New();
                        int npoin = (geo_order + 1) * (geo_order + 2) * (geo_order + 3) / 6;
                        tetr->GetPointIds()->SetNumberOfIds(npoin);
                        tetr->GetPoints()->SetNumberOfPoints(npoin);
                        tetr->Initialize();
                        reference_simplices[geo_order] = tetr;
                    }
                    break;
            }
        }

        // setup unstructured grid
        auto vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
#ifdef ICEICLE_USE_MPI
        // set VTK to use our MPI communicator
        vtkSmartPointer<vtkMPICommunicator> vtk_comm = vtkSmartPointer<vtkMPICommunicator>::New();
        MPI_Comm mpi_comm = MPI_COMM_WORLD;
        vtkMPICommunicatorOpaqueComm vtk_opaque_comm(&mpi_comm);
        vtk_comm->InitializeExternal(&vtk_opaque_comm);
        vtkSmartPointer<vtkMPIController> vtk_mpi_ctrl = vtkSmartPointer<vtkMPIController>::New();
        vtk_mpi_ctrl->SetCommunicator(vtk_comm);
        writer->SetController(vtk_mpi_ctrl);
#endif

        vtk_grid->Allocate(mesh.nelem());

        // ===================
        // = Write the Nodes =
        // ===================
        // This is done by looping through the elements, getting the VTK reference cell 
        // then interpolating to vtk's reference domain points

        IDX ivtk_pt = 0;
        vtkNew<vtkPoints> vtk_coord;
        for(IDX ielem = 0; ielem < mesh.nelem(); ++ielem){
            switch(mesh.el_transformations[ielem]->domain_type){
                case DOMAIN_TYPE::HYPERCUBE:
                    auto ref_cell = reference_hypercubes[mesh.el_transformations[ielem]->order];
                    double* coords = ref_cell->GetParametricCoords();
                    std::vector<vtkIdType> pts{};
                    for(int ipoin = 0; ipoin < ref_cell->GetNumberOfPoints(); ++ipoin){
                        double* coords_start = coords + 3 * ipoin; // always 3d 
                        MATH::GEOMETRY::Point<T, ndim> refpt{};
                        for(int idim = 0; idim < ndim; ++idim){
                            // shift from [0, 1] to [-1, 1]
                            refpt[idim] = coords_start[idim] * 2.0 - 1.0;
                        }
                        std::vector<MATH::GEOMETRY::Point<T, ndim>> el_coord(
                                mesh.coord_els.rowsize(ielem));
                        std::ranges::copy(mesh.coord_els.rowspan(ielem), el_coord.begin());
                        auto act_pt = mesh.el_transformations[ielem]->transform(
                            std::span{el_coord}, refpt);
                        double vtk_point[3];
                        for(int idim = 0; idim < std::min(3, ndim); ++idim)
                            vtk_point[idim] = act_pt[idim];
                        for(int idim = ndim; idim < 3; ++idim)
                            vtk_point[idim] = 0.0;
                        vtk_coord->InsertPoint(ivtk_pt, vtk_point);
                        pts.push_back(ivtk_pt);
                        ++ivtk_pt;
                    }
                    vtk_grid->InsertNextCell(ref_cell->GetCellType(),
                            ref_cell->GetNumberOfPoints(), pts.data());
                break;
            }
        }
        vtk_grid->SetPoints(vtk_coord);

        // Create an Array to represent the MPI rank for each cell 
        vtkSmartPointer<vtkIntArray> rank_array = vtkSmartPointer<vtkIntArray>::New();
        rank_array->SetName("MPI rank");
        rank_array->SetNumberOfComponents(1);
        rank_array->SetNumberOfTuples(vtk_grid->GetCells()->GetNumberOfCells());
        for(IDX ielem = 0; ielem < mesh.nelem(); ++ielem) {
            rank_array->SetValue(ielem, mpi::mpi_world_rank());
        }
        vtk_grid->GetCellData()->AddArray(rank_array);

        writer->SetNumberOfPieces(mpi::mpi_world_size());
        writer->SetStartPiece(mpi::mpi_world_rank());
        writer->SetEndPiece(mpi::mpi_world_rank());

        writer->SetInputData(vtk_grid);
        writer->SetFileName("test.pvtu");
        writer->SetDataModeToAscii();
        writer->Write();
    }
}
