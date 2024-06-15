#include "gtest/gtest.h"
#include "iceicle/transformations/polytope_transformations.hpp"

using namespace iceicle;
using namespace polytope;

//  === Element Topologies ===
static constexpr bitset<1> segment_t{"0"};
static constexpr tcode<1> segmentb_t{"1"};

static constexpr tcode<2> tri_a_t("00");
static constexpr tcode<2> tri_b_t("01");
static constexpr tcode<2> quad_a_t("10");
static constexpr tcode<2> quad_b_t("11");

static constexpr tcode<3> tet_a_t{"000"};
static constexpr tcode<3> tet_b_t{"001"};
static constexpr tcode<3> pyra_a_t{"010"};
static constexpr tcode<3> pyra_b_t{"011"};
static constexpr tcode<3> prism_a_t{"100"};
static constexpr tcode<3> prism_b_t{"101"};
static constexpr tcode<3> hexa_a_t{"110"};
static constexpr tcode<3> hexa_b_t{"111"};

TEST(test_polytope, test_vertices){

    // === 1D elements ===

    ASSERT_EQ(1, get_ndim(segment_t));
    ASSERT_EQ(2, n_vert(segment_t));
    ASSERT_EQ(2, n_vert(segmentb_t));

    {
        auto expected_vlist = std::array{
            bitset<1>{"0"},
            bitset<1>{"1"}
        };
        ASSERT_EQ(gen_vert<segment_t>(), expected_vlist );
        ASSERT_EQ(gen_vert<segmentb_t>(), expected_vlist );
    }

    // ==== 2D elements ===
    { // triangles
        auto expected_vlist = std::array{
            vcode<2>{"00"},
            vcode<2>{"01"},
            vcode<2>{"10"}
        };
        ASSERT_EQ(gen_vert<tri_a_t>(), expected_vlist );
        ASSERT_EQ(gen_vert<tri_b_t>(), expected_vlist );
    }

    { // quads
        auto expected_vlist = std::array{
            vcode<2>{"00"},
            vcode<2>{"01"},
            vcode<2>{"10"},
            vcode<2>{"11"}
        };
        ASSERT_EQ(gen_vert<quad_a_t>(), expected_vlist );
        ASSERT_EQ(gen_vert<quad_b_t>(), expected_vlist );
    }

    // === 3D elements ===
    { // tetrahedron
        auto expected_vlist = std::array{
            vcode<3>{"000"},
            vcode<3>{"001"},
            vcode<3>{"010"},
            vcode<3>{"100"}
        };
        ASSERT_EQ(gen_vert<tet_a_t>(), expected_vlist );
        ASSERT_EQ(gen_vert<tet_b_t>(), expected_vlist );
    }

    { // pyramid 
        auto expected_vlist = std::array{
            vcode<3>{"000"},
            vcode<3>{"001"},
            vcode<3>{"010"},
            vcode<3>{"011"},
            vcode<3>{"100"}
        };
        ASSERT_EQ(gen_vert<pyra_a_t>(), expected_vlist );
        ASSERT_EQ(gen_vert<pyra_b_t>(), expected_vlist );
    }

    { // triangle prism
        auto expected_vlist = std::array{
            vcode<3>{"000"},
            vcode<3>{"001"},
            vcode<3>{"010"},
            vcode<3>{"100"},
            vcode<3>{"101"},
            vcode<3>{"110"}
        };
        ASSERT_EQ(gen_vert<prism_a_t>(), expected_vlist );
        ASSERT_EQ(gen_vert<prism_b_t>(), expected_vlist );
    }

    { // hexahedron
        auto expected_vlist = std::array{
            vcode<3>{"000"},
            vcode<3>{"001"},
            vcode<3>{"010"},
            vcode<3>{"011"},
            vcode<3>{"100"},
            vcode<3>{"101"},
            vcode<3>{"110"},
            vcode<3>{"111"}
        };
        ASSERT_EQ(gen_vert<hexa_a_t>(), expected_vlist );
        ASSERT_EQ(gen_vert<hexa_b_t>(), expected_vlist );
    }
}
