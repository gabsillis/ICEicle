/** 
 * @file segment.hpp 
 * @author Gianni Absillis (gabsill@ncsu.edu)
 * @brief A line segment
 * 
 */
#pragma once
#include "iceicle/fe_definitions.hpp"
#include <iceicle/geometry/geo_element.hpp>
#include <iceicle/transformations/SegmentTransformation.hpp>

namespace iceicle {

    template<typename T, typename IDX>
    class Segment final : public GeometricElement<T, IDX, 1> {
        private:
        static constexpr int ndim = 1;
        static constexpr int nnodes = 2;

        using Point = MATH::GEOMETRY::Point<T, ndim>;
        using HessianType = NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim, ndim>;

        IDX node_idxs[nnodes];

        public:

        static inline transformations::SegmentTransformation<T, IDX> transformation{};

        Segment(IDX node1, IDX node2) {
            node_idxs[0] = node1;
            node_idxs[1] = node2;
        }

        constexpr int n_nodes() const override { return nnodes; }

        constexpr DOMAIN_TYPE domain_type() const noexcept override { return DOMAIN_TYPE::HYPERCUBE; }

        constexpr int geometry_order() const noexcept override { return 1; }

        const IDX *nodes() const override { return node_idxs; }

        void transform(NodeArray<T, ndim> &node_coords, const Point &pt_ref, Point &pt_phys) const  override {
            return transformation.transform(node_coords, node_idxs, pt_ref, pt_phys);        
        }
        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim> Jacobian(
            NodeArray<T, ndim> &node_coords,
            const Point &xi
        ) const override 
        { return transformation.Jacobian(node_coords, node_idxs, xi); }


        HessianType Hessian(
            NodeArray<T, ndim> &node_coords,
            const Point &xi
        ) const override {
            HessianType ret{{{{0}}}};;
            return ret;
        }

        auto regularize_interior_nodes(
            NodeArray<T, ndim>& coord /// [in/out] the node coordinates array 
        ) const -> void override {
            // do nothing: no interior nodes
        }
    };
}
