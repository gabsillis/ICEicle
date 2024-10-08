/**
 * @author Gianni Absillis (gabsill@ncsu.edu)
 * @brief cglayout represents the memoery layout of
 * a CG representation of a vector-valued fe_function
 */
#pragma once
#include "iceicle/mesh/mesh.hpp"
#include <iceicle/fe_function/layout_enums.hpp>

#include <type_traits>

namespace iceicle {

    namespace impl {

        template<class T, class IDX, int ndim>
        static std::optional<AbstractMesh<T, IDX, ndim>> empty_mesh_opt = std::nullopt;

        template<class T, class IDX, int ndim>
        constexpr
        auto get_empty_mesh() -> AbstractMesh<T, IDX, ndim>& {
            if(empty_mesh_opt<T, IDX, ndim>){
                return empty_mesh_opt<T, IDX, ndim>.value();
            } else {
                empty_mesh_opt<T, IDX, ndim> = std::optional{AbstractMesh<T, IDX, ndim>{}};
                return empty_mesh_opt<T, IDX, ndim>.value();
            }
        }
    }

    /// @brief map degrees of freedom (dofs) for isoparametric CG spaces
    /// This represents a map from the pair (ielem, ildof) -> to a global dof 
    /// where ielem is the element index and ildof is the local dof on the element
    template< class T, class IDX, int ndim >
    class cg_dof_map {
    public:

        // ============
        // = Typedefs =
        // ============

        using index_type = IDX;
        using size_type = std::make_unsigned_t<index_type>;


        // ================
        // = Data Members =
        // ================

        /// reference to the mesh that we map an isoparametric space to
        AbstractMesh<T, IDX, ndim>& mesh;

        // ================
        // = Constructors =
        // ================

        constexpr 
        cg_dof_map() : mesh{impl::get_empty_mesh<T, IDX, ndim>()} {}

        constexpr 
        cg_dof_map(AbstractMesh<T, IDX, ndim>& mesh) : mesh{mesh} {}

        // === Default nothrow copy and nothrow move semantics ===
        constexpr cg_dof_map(const cg_dof_map<T, IDX, ndim>& other) noexcept 
        : mesh{other.mesh} {}
        constexpr cg_dof_map(cg_dof_map<T, IDX, ndim>&& other) noexcept 
        : mesh{other.mesh} {}

        constexpr cg_dof_map& operator=(const cg_dof_map<T, IDX, ndim>& other) noexcept {
            mesh = other.mesh;
            return *this;
        };
        constexpr cg_dof_map& operator=(cg_dof_map<T, IDX, ndim>&& other) noexcept {
            mesh = other.mesh;
            return *this;
        }

        // =============
        // = Accessors =
        // =============

        /** 
         * @brief Convert element index and local degree of freedom index 
         * to the global degree of freedom index 
         * @param ielem the element index 
         * @param idof the local degree of freedom index
         */
        constexpr index_type operator[](index_type ielem, index_type idof) 
            const noexcept { return mesh.elements[ielem]->nodes(idof); }

        // ===========
        // = Utility =
        // ===========

        /** @brief get the size requirement for all degrees of freedom given
         * the number of vector components per dof 
         * @param nv_comp th number of vector components per dof 
         * @return the size requirement
         */
        constexpr size_type calculate_size_requirement( index_type nv_comp ) const noexcept {
            return size() * nv_comp;
        }

        /**
         * @brief calculate the largest size requirement for a single element 
         * @param nv_comp the number of vector components per dof 
         * @return the maximum size requirement 
         */
        constexpr size_type max_el_size_reqirement( index_type nv_comp ) const noexcept {
            size_type max = 0;
            for(auto& element_ptr : mesh.elements){
                max = std::max(element_ptr->n_nodes() * nv_comp, max);
            }
            return max;
        }

        /**
         * @brief get the number of degrees of freedom at the given element index 
         * @param elidx the index of the element to get the ndofs for 
         * @return the number of degrees of freedom 
         */
        [[nodiscard]] constexpr size_type ndof_el( index_type elidx ) const noexcept {
            return mesh.elements[elidx]->n_nodes();
        }

        /** @brief get the number of elements represented in the map */
        [[nodiscard]] constexpr size_type nelem() const noexcept { return mesh.elements.size(); }

        /** @brief get the size of the global degree of freedom index space represented by this map */
        [[nodiscard]] constexpr size_type size() const noexcept { return static_cast<size_type>(mesh.nodes.size()); }
    };

    // Deduction guide
    template< class T, class IDX, int ndim >
    cg_dof_map(const AbstractMesh<T, IDX, ndim>&) -> cg_dof_map<T, IDX, ndim>;

}
