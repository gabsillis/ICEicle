#include "iceicle/disc/SpacetimeConnection.hpp"
#include "iceicle/geometry/face.hpp"
#include "iceicle/mesh/mesh.hpp"
#include "gtest/gtest.h"

using namespace iceicle;

TEST(test_spacetime_utils, test_node_connectivity){
    using T = double;
    using IDX = int;
    static constexpr int ndim = 2;
    AbstractMesh<T, IDX, ndim> mesh {
        {0.0, 0.0},
        {1.0, 1.0},
        {4, 4}, 
        1,
        {
            BOUNDARY_CONDITIONS::DIRICHLET,
            BOUNDARY_CONDITIONS::SPACETIME_PAST,
            BOUNDARY_CONDITIONS::DIRICHLET,
            BOUNDARY_CONDITIONS::SPACETIME_FUTURE
        }
    };

    std::map<IDX, IDX> node_connectivity{compute_st_node_connectivity(mesh, mesh)};
    ASSERT_EQ(node_connectivity[0], 20);
    ASSERT_EQ(node_connectivity[1], 21);
    ASSERT_EQ(node_connectivity[2], 22);
    ASSERT_EQ(node_connectivity[3], 23);
    ASSERT_EQ(node_connectivity[4], 24);

    
}
