/**
 * @author Gianni Absillis (gabsill@ncsu.edu)
 * @brief common definitions for memory layouts
 */
#pragma once
#include <cstdlib>
#include <span>

namespace iceicle {

    /**
     * @brief defines how degrees of freedom are organized wrt vector components
     */
    enum LAYOUT_VECTOR_ORDER {
        /**
         * degrees of freedom are the left most index in row major order setting 
         * so each degree of freedom will have a stride of the number of vector components 
         */
        DOF_LEFT = 0,
        /** degrees of freedom are the rightmost index in row major order setting
         * indices are dof fastest and in large chunks for each vector component 
         */
        DOF_RIGHT
    };

    /** @brief tag to specify that the number of vector component is a runtime parameter */
    static constexpr std::size_t dynamic_ncomp = std::dynamic_extent;

    /**
     * @brief determine if a size field represents a dynamic extent
     */
    template<std::size_t ncomp>
    struct is_dynamic_size {

        /// the boolean value
        inline static constexpr bool value = (ncomp == std::dynamic_extent);
    };


    /**
     * @brief index into a fespan 
     * collects the 3 required indices 
     * all indices default to zero
     */
    struct fe_index {
        /// the element index 
        std::size_t iel = 0;

        /// the local degree of freedom index 
        std::size_t idof = 0;

        /// the index of the vector component
        std::size_t iv = 0;
    };

    /**
     * @brief index into a element local elspan
     * collects the 2 required indices 
     * all indices default to zero 
     */
    struct compact_index {
        /// the local degree of freedom index 
        std::size_t idof = 0;

        /// the index of the vector component 
        std::size_t iv = 0;
    };

    /**
     * @brief the extents of the index space for 
     * the mutidimensional compact_index 
     */
    template<std::size_t ncomp>
    struct compact_index_extents {
        std::size_t ndof = 0;

        std::size_t nv = (is_dynamic_size<ncomp>::value) ?
            0 : ncomp;

        /// get the ncomp template argument
        static constexpr std::size_t get_ncomp() noexcept { return ncomp; }
    };

    /**
     * @brief index into a global nodal structure
     * collects the two required indices 
     * indices default to zero 
     */
    struct gnode_index {
        /// the node index 
        std::size_t idof = 0;

        /// the index of the vector component
        std::size_t iv = 0;
    };

    /**
     * @brief struct concept to tell if the data for an element can be 
     * copied out of the global fespan in a contiguous block
     */
    template< class local_layout, class global_layout >
    struct is_equivalent_el_layout {
        static constexpr bool value = false;
    };

    /**
     * @brief struct concept to tell if the data in a dofspan 
     * can be copied out to the global fespan in a contiguous block 
     */
    template< class local_span, class global_span>
    struct has_equivalent_el_layout 
    : is_equivalent_el_layout<
        typename local_span::layout_type,
        typename global_span::layout_type
    >{};

}
