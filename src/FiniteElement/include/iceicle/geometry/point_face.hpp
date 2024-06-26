/**
 * @file point_face.hpp
 * @author Gianni Absillis (gabsill@ncsu.edu)
 * @brief A 1D face type
 *
 */
#pragma once
#include "iceicle/fe_definitions.hpp"
#include <iceicle/geometry/face.hpp>

namespace iceicle {
    
    template<typename T, typename IDX>
    class PointFace final : public Face<T, IDX, 1> {
        static constexpr int ndim = 1;
        using Point = MATH::GEOMETRY::Point<T, ndim>;
        using FacePoint = MATH::GEOMETRY::Point<T, ndim - 1>;

        // === Index data ===
        /// the node corresponding to this face
        IDX node;

        // === Precomputation ===
        /// the normal vector
        T normal;
        /// the face area
        T area = 1.0;

        public:
        // === Constructors ===
        PointFace(
            IDX elemL,
            IDX elemR,
            int faceNrL,
            int faceNrR,
            IDX node,
            bool positiveNormal,
            BOUNDARY_CONDITIONS bctype = BOUNDARY_CONDITIONS::INTERIOR,
            int bcflag = 0
        ) : Face<T, IDX, 1>(
                elemL, elemR,
                faceNrL * FACE_INFO_MOD, 
                faceNrR * FACE_INFO_MOD,
                bctype, bcflag
            ),
            node(node)
        {
           if(positiveNormal) normal = 1.0;
           else normal = -1.0;
        }

        constexpr DOMAIN_TYPE domain_type() const override { return DOMAIN_TYPE::HYPERCUBE; }

        constexpr int geometry_order() const noexcept override { return 1; }

        T getArea() const override { return area; }

        void transform(
            const FacePoint &s,
            NodeArray<T,ndim> &nodeCoords,
            Point& result
        ) const override {
            result[0] = nodeCoords[node][0];
        }

        void transform_xiL(
            const FacePoint &s,
            Point& result
        ) const override {
            result[0] = -1.0;
        }

        void transform_xiR(
            const FacePoint &s,
            Point& result
        ) const override {
            result[0] = 1.0;
        }


        NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim-1> Jacobian(
            NodeArray<T, ndim> &node_coords,
            const FacePoint &s
        ) const override {
            NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim-1> ret;
            ret[0][0] = normal; // use extra space defined in the Tensor for zero size
            return ret;
        }


        virtual
        T rootRiemannMetric(
            NUMTOOL::TENSOR::FIXED_SIZE::Tensor<T, ndim, ndim - 1> &J,
            const FacePoint &s
        ) const override { return 1.0; }

        int n_nodes() const override { return 1; }

        const IDX *nodes() const override { return &node; }
    };
}
