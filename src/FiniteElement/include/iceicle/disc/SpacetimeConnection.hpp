
#pragma once

#include "iceicle/element/TraceSpace.hpp"
#include "iceicle/fe_function/fespan.hpp"
#include "iceicle/fespace/fespace.hpp"
#include "iceicle/geometry/face.hpp"
#include "iceicle/mesh/mesh.hpp"
#include <map>
#include <list>

namespace MESH {

    /**
     *
     *  NOTE: we assume the time dimension is always the last
     */
    template<class T, class IDX, int ndim>
    auto compute_st_node_connectivity(
        MESH::AbstractMesh<T, IDX, ndim> &mesh_past,
        MESH::AbstractMesh<T, IDX, ndim> &mesh_current
    ) -> std::map<IDX, IDX> {
        static_assert(ndim > 1, "Assumes at least 2D spacetime mesh");
        using Face = ELEMENT::Face<T, IDX, ndim>;
        std::vector<bool> past_nodes_connected(mesh_past.nodes.n_nodes(), false);
        std::vector<bool> current_nodes_connected(mesh_current.nodes.n_nodes(), false);

        // look for nodes in SPACETIME_FUTURE boundary in the past mesh
        for(IDX ibface = mesh_past.bdyFaceStart; ibface < mesh_past.bdyFaceEnd; ++ibface){
            const Face& face = *(mesh_past.faces[ibface]);
            if(face.bctype == ELEMENT::BOUNDARY_CONDITIONS::SPACETIME_FUTURE){
                for(IDX inode : face.nodes_span()){
                    past_nodes_connected[inode] = true;
                }
            }
        }

        // look for nodes in SPACETIME_PAST boundary in the current mesh
        for(IDX ibface = mesh_current.bdyFaceStart; ibface < mesh_current.bdyFaceEnd; ++ibface){
            const Face& face = *(mesh_current.faces[ibface]);
            if(face.bctype == ELEMENT::BOUNDARY_CONDITIONS::SPACETIME_PAST){
                for(IDX inode : face.nodes_span()){
                    current_nodes_connected[inode] = true;
                }
            }
        }

        // create array of the node indices on the past mesh that are on the ST boundary
        std::vector<IDX> past_nodes{};
        for(IDX inode = 0; inode < mesh_past.nodes.n_nodes(); ++inode){
            if(past_nodes_connected[inode]) past_nodes.push_back(inode);
        }

        // find connected nodes 
        std::map<IDX, IDX> curr_to_past_nodes;
        for(IDX inode_curr = 0; inode_curr < mesh_current.nodes.n_nodes(); ++inode_curr){
            if(current_nodes_connected[inode_curr]){
                // search for a node in the past_nodes list that is close enough in all coordinates except the last 
                // NOTE: we assume the time dimension is always the last
                for(IDX inode_past : past_nodes){
                    bool all_same = true;
                    for(int idim = 0; idim < ndim - 1; ++idim){
                        if(std::abs(mesh_current.nodes[inode_curr][idim] - mesh_past.nodes[inode_past][idim]) > 1e-8){
                            all_same = false;
                            break;
                        }
                    }

                    if(all_same){
                        curr_to_past_nodes[inode_curr] = inode_past;
                        break;
                    }
                }
            }
        }

        return curr_to_past_nodes;
    }
}

namespace DISC {
    template<
        class T,
        class IDX,
        int ndim,
        class upastLayout
    >
    class SpacetimeConnection {

        FE::FESpace<T, IDX, ndim> &fespace_past;
        FE::FESpace<T, IDX, ndim> &fespace_current;
        FE::fespan<T, upastLayout> u_past;

        /// map the index of nodes in the current fespace 
        /// to nodes in the past fespace
        std::map<IDX, IDX> &curr_to_past_nodes;

        public:

        SpacetimeConnection(
            FE::FESpace<T, IDX, ndim> &fespace_past,
            FE::FESpace<T, IDX, ndim> &fespace_current,
            FE::fespan<T, upastLayout> &u_past,
            std::map<IDX, IDX> &curr_to_past_nodes
        ) noexcept : fespace_past{fespace_past}, fespace_current{fespace_current}, 
                     u_past{u_past}, curr_to_past_nodes{curr_to_past_nodes}
        {
            using TraceSpace = ELEMENT::TraceSpace<T, IDX, ndim>;

            // === find attached traces ===
            
            // boundary face indexes that still need to be connected
            std::list<IDX> past_bface_idxs{};
            for(int ibface = fespace_past.bdy_trace_start; ibface < fespace_past.bdy_trace_end; ++ibface) {
                TraceSpace past_trace = fespace_past.traces[ibface];

                // from the perspective of past fespace, the current fespace faces are connected at SPACETIME_FUTURE
                if(past_trace.face.bctype == ELEMENT::BOUNDARY_CONDITIONS::SPACETIME_FUTURE) {
                    past_bface_idxs.push_back(ibface);
                }
            }

            for(TraceSpace &current_trace : fespace_current.get_boundary_traces()) {
                
            }

        }
    };

}
