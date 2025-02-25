#include "gtest/gtest.h"
#include "iceicle/basis/lagrange_1d.hpp"
#include "iceicle/transformations/polytope_transformations.hpp"
#include <Numtool/fixed_size_tensor.hpp>
#include <fmt/ranges.h>

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

TEST(test_polytope, test_extrusion_parities) {
    {
        ecode<3> e{"011"};
        vcode<3> v{"100"};
        ASSERT_TRUE(extrusion_parity(e, v));
    }

    {
        ecode<3> e{"011"};
        vcode<3> v{"110"};
        ASSERT_FALSE(extrusion_parity(e, v));
    }

    // test the notion of ccw = outward normal
    { // consider: x = 0 face of unit cube

        // if we choose v = 000 and extrude y then z
        // right hand rule tells us the normal points into the domain
        ecode<3> e{"110"};
        vcode<3> v1{"000"};

        ASSERT_NE(
            extrusion_parity(~e, v1), // dual extrusion 
            hodge_extrusion_parity(e, v1)
        );

        // if we choose v = 010 or v = 100 and extrude y then z
        // right hand rule tells use the normal points out of the domain (ccw -> outward normal)
        vcode<3> v2{"010"};
        vcode<3> v3{"100"};
        ASSERT_EQ(
            extrusion_parity(~e, v2), // dual extrusion 
            hodge_extrusion_parity(e, v2)
        );
        ASSERT_EQ(
            extrusion_parity(~e, v3), // dual extrusion 
            hodge_extrusion_parity(e, v3)
        );

    }
}

// TEST(test_polytope, test_extrusion_vertices){
//     {
//         constexpr tcode<2> tria_t{"00"};
//         {
//   constexpr ecode<2> e{"10"};
//             constexpr vcode<2> v{"01"};
//             ASSERT_EQ(
//                 (facet_vertices<tria_t, e, v>()),
//                 (std::array{vcode<2>{"01"}, vcode<2>{"10"}})
//             );
//         }
//         {
//             constexpr ecode<2> e{"01"};
//             constexpr vcode<2> v{"10"};
//             ASSERT_EQ(
//                 (facet_vertices<tria_t, e, v>()),
//                 (std::array{vcode<2>{"10"}, vcode<2>{"01"}})
//             );
//         }
//     }
//     {
//         constexpr tcode<3> pyra_t{"011"};
//         
//         {
//             constexpr ecode<3> e{"000"};
//             constexpr vcode<3> v{"010"};
//             ASSERT_EQ(
//                 (facet_vertices<pyra_t, e, v>()),
//                 (std::array{vcode<3>{"010"}})
//             );
//             ASSERT_EQ((extrusion_topology<pyra_t, e>()), (tcode<0>{}));
//         }
// 
//         {
//             constexpr ecode<3> e{"110"};
//             constexpr vcode<3> v{"011"};
//             ASSERT_EQ(
//                 (facet_vertices<pyra_t, e, v>()),
//                 (std::array{vcode<3>{"011"}, vcode<3>{"001"}, vcode<3>{"100"}})
//             );
//             ASSERT_EQ((extrusion_topology<pyra_t, e>()), (tcode<2>{"01"}));
//         }
//     }
//     {
//         constexpr tcode<3> pyra_t{"010"};
//         
//         {
//             constexpr ecode<3> e{"000"};
//             constexpr vcode<3> v{"010"};
//             ASSERT_EQ(
//                 (facet_vertices<pyra_t, e, v>()),
//                 (std::array{vcode<3>{"010"}})
//             );
//         }
// 
//         {
//             constexpr ecode<3> e{"110"};
//             constexpr vcode<3> v{"011"};
//             ASSERT_EQ(
//                 (facet_vertices<pyra_t, e, v>()),
//                 (std::array{vcode<3>{"011"}, vcode<3>{"001"}, vcode<3>{"100"}})
//             );
//             ASSERT_EQ((extrusion_topology<pyra_t, e>()), (tcode<2>{"00"}));
//         }
//         {
//             constexpr ecode<3> e{"101"};
//             constexpr vcode<3> v{"001"};
//             auto vertices = facet_vertices<pyra_t, e, v>();
//             for(auto vert : vertices){
//                 std::cout << vert.to_string() << " " << std::endl;
//             }
//         }
//     }
//     {
//         constexpr tcode<3> tri_prism_t{"101"};
//         {
//             constexpr ecode<3> e{"101"};
//             constexpr vcode<3> v{"011"};
//             ASSERT_EQ(
//                 (facet_vertices<tri_prism_t, e, v>()),
//                 (std::array{vcode<3>{"011"}, vcode<3>{"101"}, vcode<3>{"010"}, vcode<3>{"001"}})
//             );
//             
//         }
//     }
// }
//
TEST(test_polytope, test_n_facets){
    // Vertices
    ASSERT_EQ(n_facets<0>(segment_t), 2);
    ASSERT_EQ(n_facets<0>(segmentb_t), 2);
    ASSERT_EQ(n_facets<0>(tri_a_t), 3);
    ASSERT_EQ(n_facets<0>(tri_b_t), 3);
    ASSERT_EQ(n_facets<0>(quad_a_t), 4);
    ASSERT_EQ(n_facets<0>(quad_b_t), 4);
    ASSERT_EQ(n_facets<0>(tet_a_t), 4);
    ASSERT_EQ(n_facets<0>(tet_b_t), 4);
    ASSERT_EQ(n_facets<0>(prism_a_t), 6);
    ASSERT_EQ(n_facets<0>(prism_b_t), 6);
    ASSERT_EQ(n_facets<0>(pyra_a_t), 5);
    ASSERT_EQ(n_facets<0>(pyra_b_t), 5);
    ASSERT_EQ(n_facets<0>(hexa_a_t), 8);
    ASSERT_EQ(n_facets<0>(hexa_b_t), 8);

    // lines
    ASSERT_EQ(n_facets<1>(segment_t), 1);
    ASSERT_EQ(n_facets<1>(segmentb_t), 1);
    ASSERT_EQ(n_facets<1>(tri_a_t), 3);
    ASSERT_EQ(n_facets<1>(tri_b_t), 3);
    ASSERT_EQ(n_facets<1>(quad_a_t), 4);
    ASSERT_EQ(n_facets<1>(quad_b_t), 4);
    ASSERT_EQ(n_facets<1>(tet_a_t), 6);
    ASSERT_EQ(n_facets<1>(tet_b_t), 6);
    ASSERT_EQ(n_facets<1>(prism_a_t), 9);
    ASSERT_EQ(n_facets<1>(prism_b_t), 9);
    ASSERT_EQ(n_facets<1>(pyra_a_t), 8);
    ASSERT_EQ(n_facets<1>(pyra_b_t), 8);
    ASSERT_EQ(n_facets<1>(hexa_a_t), 12);
    ASSERT_EQ(n_facets<1>(hexa_b_t), 12);

    // surfaces
    ASSERT_EQ(n_facets<2>(segment_t), 0);
    ASSERT_EQ(n_facets<2>(segmentb_t), 0);
    ASSERT_EQ(n_facets<2>(tri_a_t), 1);
    ASSERT_EQ(n_facets<2>(tri_b_t), 1);
    ASSERT_EQ(n_facets<2>(quad_a_t), 1);
    ASSERT_EQ(n_facets<2>(quad_b_t), 1);
    ASSERT_EQ(n_facets<2>(tet_a_t), 4);
    ASSERT_EQ(n_facets<2>(tet_b_t), 4);
    ASSERT_EQ(n_facets<2>(prism_a_t), 5);
    ASSERT_EQ(n_facets<2>(prism_b_t), 5);
    ASSERT_EQ(n_facets<2>(pyra_a_t), 5);
    ASSERT_EQ(n_facets<2>(pyra_b_t), 5);
    ASSERT_EQ(n_facets<2>(hexa_a_t), 6);
    ASSERT_EQ(n_facets<2>(hexa_b_t), 6);

    // volumes
    ASSERT_EQ(n_facets<3>(segment_t), 0);
    ASSERT_EQ(n_facets<3>(segmentb_t), 0);
    ASSERT_EQ(n_facets<3>(tri_a_t), 0);
    ASSERT_EQ(n_facets<3>(tri_b_t), 0);
    ASSERT_EQ(n_facets<3>(quad_a_t), 0);
    ASSERT_EQ(n_facets<3>(quad_b_t), 0);
    ASSERT_EQ(n_facets<3>(tet_a_t), 1);
    ASSERT_EQ(n_facets<3>(tet_b_t), 1);
    ASSERT_EQ(n_facets<3>(prism_a_t), 1);
    ASSERT_EQ(n_facets<3>(prism_b_t), 1);
    ASSERT_EQ(n_facets<3>(pyra_a_t), 1);
    ASSERT_EQ(n_facets<3>(pyra_b_t), 1);
    ASSERT_EQ(n_facets<3>(hexa_a_t), 1);
    ASSERT_EQ(n_facets<3>(hexa_b_t), 1);
}

TEST(test_polytope, test_facet_definitions){

    // ===========
    // = Segment =
    // ===========
    {
        // vertices
        std::vector<facet<0>> vertices = get_facets<0>(segment_t);
        ASSERT_EQ(vertices.size(), 2);
        ASSERT_EQ(vertices[0].t, tcode<0>{});
        ASSERT_EQ(vertices[0].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[0].vertex_indices[0], 0);
        ASSERT_EQ(vertices[1].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[1].vertex_indices[0], 1);

        // lines
        std::vector<facet<1>> lines = get_facets<1>(segment_t);
        ASSERT_EQ(lines.size(), 1);
        ASSERT_EQ(lines[0].t, tcode<1>{"0"});
        ASSERT_EQ(lines[0].vertex_indices, (std::vector<std::size_t>{0, 1}));

        // surfaces
        std::vector<facet<2>> surfaces = get_facets<2>(segment_t);
        ASSERT_EQ(surfaces.size(), 0);

        // volumes 
        std::vector<facet<3>> volumes = get_facets<3>(segment_t);
        ASSERT_EQ(volumes.size(), 0);
    }
    {
        // vertices
        std::vector<facet<0>> vertices = get_facets<0>(segmentb_t);
        ASSERT_EQ(vertices.size(), 2);
        ASSERT_EQ(vertices[0].t, tcode<0>{});
        ASSERT_EQ(vertices[0].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[0].vertex_indices[0], 0);
        ASSERT_EQ(vertices[1].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[1].vertex_indices[0], 1);

        // lines
        std::vector<facet<1>> lines = get_facets<1>(segmentb_t);
        ASSERT_EQ(lines.size(), 1);
        ASSERT_EQ(lines[0].t, tcode<1>{"1"});
        ASSERT_EQ(lines[0].vertex_indices, (std::vector<std::size_t>{0, 1}));

        // surfaces
        std::vector<facet<2>> surfaces = get_facets<2>(segmentb_t);
        ASSERT_EQ(surfaces.size(), 0);

        // volumes 
        std::vector<facet<3>> volumes = get_facets<3>(segmentb_t);
        ASSERT_EQ(volumes.size(), 0);
    }

    // ============
    // = Triangle =
    // ============
    {
        // vertices
        std::vector<facet<0>> vertices = get_facets<0>(tri_a_t);
        ASSERT_EQ(vertices.size(), 3);
        ASSERT_EQ(vertices[0].t, tcode<0>{});
        ASSERT_EQ(vertices[0].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[0].vertex_indices, std::vector<std::size_t>{0});
        ASSERT_EQ(vertices[1].t, tcode<0>{});
        ASSERT_EQ(vertices[1].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[1].vertex_indices, std::vector<std::size_t>{1});
        ASSERT_EQ(vertices[2].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[2].vertex_indices, std::vector<std::size_t>{2});

        // lines
        std::vector<facet<1>> lines = get_facets<1>(tri_a_t);
        ASSERT_EQ(lines.size(), 3);
        ASSERT_EQ(lines[0].t, tcode<1>{"0"});
        ASSERT_EQ(lines[0].vertex_indices, (std::vector<std::size_t>{0, 1}));
        ASSERT_EQ(lines[1].t, tcode<1>{"0"});
        ASSERT_EQ(lines[1].vertex_indices, (std::vector<std::size_t>{0, 2}));
        ASSERT_EQ(lines[2].t, tcode<1>{"0"});
        ASSERT_EQ(lines[2].vertex_indices, (std::vector<std::size_t>{1, 2}));

        // surfaces
        std::vector<facet<2>> surfaces = get_facets<2>(tri_a_t);
        ASSERT_EQ(surfaces.size(), 1);
        ASSERT_EQ(surfaces[0].t, tcode<2>{"00"});
        ASSERT_EQ(surfaces[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2}));

        // volumes 
        std::vector<facet<3>> volumes = get_facets<3>(tri_a_t);
        ASSERT_EQ(volumes.size(), 0);
    }
    {
        // vertices
        std::vector<facet<0>> vertices = get_facets<0>(tri_b_t);
        ASSERT_EQ(vertices.size(), 3);
        ASSERT_EQ(vertices[0].t, tcode<0>{});
        ASSERT_EQ(vertices[0].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[0].vertex_indices, std::vector<std::size_t>{0});
        ASSERT_EQ(vertices[1].t, tcode<0>{});
        ASSERT_EQ(vertices[1].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[1].vertex_indices, std::vector<std::size_t>{1});
        ASSERT_EQ(vertices[2].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[2].vertex_indices, std::vector<std::size_t>{2});

        // lines
        std::vector<facet<1>> lines = get_facets<1>(tri_b_t);
        ASSERT_EQ(lines.size(), 3);
        ASSERT_EQ(lines[0].t, tcode<1>{"1"});
        ASSERT_EQ(lines[0].vertex_indices, (std::vector<std::size_t>{0, 1}));
        ASSERT_EQ(lines[1].t, tcode<1>{"0"});
        ASSERT_EQ(lines[1].vertex_indices, (std::vector<std::size_t>{0, 2}));
        ASSERT_EQ(lines[2].t, tcode<1>{"0"});
        ASSERT_EQ(lines[2].vertex_indices, (std::vector<std::size_t>{1, 2}));

        // surfaces
        std::vector<facet<2>> surfaces = get_facets<2>(tri_b_t);
        ASSERT_EQ(surfaces.size(), 1);
        ASSERT_EQ(surfaces[0].t, tcode<2>{"01"});
        ASSERT_EQ(surfaces[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2}));

        // volumes 
        std::vector<facet<3>> volumes = get_facets<3>(tri_b_t);
        ASSERT_EQ(volumes.size(), 0);
    }
    // ========
    // = Quad =
    // ========
    {
        // vertices
        std::vector<facet<0>> vertices = get_facets<0>(quad_a_t);
        ASSERT_EQ(vertices.size(), 4);
        ASSERT_EQ(vertices[0].t, tcode<0>{});
        ASSERT_EQ(vertices[0].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[0].vertex_indices, std::vector<std::size_t>{0});
        ASSERT_EQ(vertices[1].t, tcode<0>{});
        ASSERT_EQ(vertices[1].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[1].vertex_indices, std::vector<std::size_t>{1});
        ASSERT_EQ(vertices[2].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[2].vertex_indices, std::vector<std::size_t>{2});
        ASSERT_EQ(vertices[3].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[3].vertex_indices, std::vector<std::size_t>{3});

        // lines
        std::vector<facet<1>> lines = get_facets<1>(quad_a_t);
        ASSERT_EQ(lines.size(), 4);
        ASSERT_EQ(lines[0].t, tcode<1>{"0"});
        ASSERT_EQ(lines[0].vertex_indices, (std::vector<std::size_t>{0, 1}));
        ASSERT_EQ(lines[1].t, tcode<1>{"0"});
        ASSERT_EQ(lines[1].vertex_indices, (std::vector<std::size_t>{2, 3}));
        ASSERT_EQ(lines[2].t, tcode<1>{"1"});
        ASSERT_EQ(lines[2].vertex_indices, (std::vector<std::size_t>{0, 2}));
        ASSERT_EQ(lines[3].t, tcode<1>{"1"});
        ASSERT_EQ(lines[3].vertex_indices, (std::vector<std::size_t>{1, 3}));

        // surfaces
        std::vector<facet<2>> surfaces = get_facets<2>(quad_a_t);
        ASSERT_EQ(surfaces.size(), 1);
        ASSERT_EQ(surfaces[0].t, tcode<2>{"10"});
        ASSERT_EQ(surfaces[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2, 3}));

        // volumes 
        std::vector<facet<3>> volumes = get_facets<3>(quad_a_t);
        ASSERT_EQ(volumes.size(), 0);
    }
    {
        // vertices
        std::vector<facet<0>> vertices = get_facets<0>(quad_b_t);
        ASSERT_EQ(vertices.size(), 4);
        ASSERT_EQ(vertices[0].t, tcode<0>{});
        ASSERT_EQ(vertices[0].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[0].vertex_indices, std::vector<std::size_t>{0});
        ASSERT_EQ(vertices[1].t, tcode<0>{});
        ASSERT_EQ(vertices[1].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[1].vertex_indices, std::vector<std::size_t>{1});
        ASSERT_EQ(vertices[2].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[2].vertex_indices, std::vector<std::size_t>{2});
        ASSERT_EQ(vertices[3].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[3].vertex_indices, std::vector<std::size_t>{3});

        // lines
        std::vector<facet<1>> lines = get_facets<1>(quad_b_t);
        ASSERT_EQ(lines.size(), 4);
        ASSERT_EQ(lines[0].t, tcode<1>{"1"});
        ASSERT_EQ(lines[0].vertex_indices, (std::vector<std::size_t>{0, 1}));
        ASSERT_EQ(lines[1].t, tcode<1>{"1"});
        ASSERT_EQ(lines[1].vertex_indices, (std::vector<std::size_t>{2, 3}));
        ASSERT_EQ(lines[2].t, tcode<1>{"1"});
        ASSERT_EQ(lines[2].vertex_indices, (std::vector<std::size_t>{0, 2}));
        ASSERT_EQ(lines[3].t, tcode<1>{"1"});
        ASSERT_EQ(lines[3].vertex_indices, (std::vector<std::size_t>{1, 3}));

        // surfaces
        std::vector<facet<2>> surfaces = get_facets<2>(quad_b_t);
        ASSERT_EQ(surfaces.size(), 1);
        ASSERT_EQ(surfaces[0].t, tcode<2>{"11"});
        ASSERT_EQ(surfaces[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2, 3}));

        // volumes 
        std::vector<facet<3>> volumes = get_facets<3>(quad_b_t);
        ASSERT_EQ(volumes.size(), 0);
    }

    // =======
    // = Tet =
    // =======
    {
        // vertices
        std::vector<facet<0>> vertices = get_facets<0>(tet_a_t);
        ASSERT_EQ(vertices.size(), 4);
        ASSERT_EQ(vertices[0].t, tcode<0>{});
        ASSERT_EQ(vertices[0].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[0].vertex_indices, std::vector<std::size_t>{0});
        ASSERT_EQ(vertices[1].t, tcode<0>{});
        ASSERT_EQ(vertices[1].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[1].vertex_indices, std::vector<std::size_t>{1});
        ASSERT_EQ(vertices[2].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[2].vertex_indices, std::vector<std::size_t>{2});
        ASSERT_EQ(vertices[3].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[3].vertex_indices, std::vector<std::size_t>{3});

        // lines
        std::vector<facet<1>> lines = get_facets<1>(tet_a_t);
        ASSERT_EQ(lines.size(), 6);
        ASSERT_EQ(lines[0].t, tcode<1>{"0"});
        ASSERT_EQ(lines[0].vertex_indices, (std::vector<std::size_t>{0, 1}));
        ASSERT_EQ(lines[1].t, tcode<1>{"0"});
        ASSERT_EQ(lines[1].vertex_indices, (std::vector<std::size_t>{0, 2}));
        ASSERT_EQ(lines[2].t, tcode<1>{"0"});
        ASSERT_EQ(lines[2].vertex_indices, (std::vector<std::size_t>{1, 2}));
        ASSERT_EQ(lines[3].t, tcode<1>{"0"});
        ASSERT_EQ(lines[3].vertex_indices, (std::vector<std::size_t>{0, 3}));
        ASSERT_EQ(lines[4].t, tcode<1>{"0"});
        ASSERT_EQ(lines[4].vertex_indices, (std::vector<std::size_t>{1, 3}));
        ASSERT_EQ(lines[5].t, tcode<1>{"0"});
        ASSERT_EQ(lines[5].vertex_indices, (std::vector<std::size_t>{2, 3}));

        // surfaces
        std::vector<facet<2>> surfaces = get_facets<2>(tet_a_t);
        ASSERT_EQ(surfaces.size(), 4);
        ASSERT_EQ(surfaces[0].t, tcode<2>{"00"});
        ASSERT_EQ(surfaces[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2}));
        ASSERT_EQ(surfaces[1].t, tcode<2>{"00"});
        ASSERT_EQ(surfaces[1].vertex_indices, (std::vector<std::size_t>{0, 1, 3}));
        ASSERT_EQ(surfaces[2].t, tcode<2>{"00"});
        ASSERT_EQ(surfaces[2].vertex_indices, (std::vector<std::size_t>{0, 2, 3}));
        ASSERT_EQ(surfaces[3].t, tcode<2>{"00"});
        ASSERT_EQ(surfaces[3].vertex_indices, (std::vector<std::size_t>{1, 2, 3}));

        // volumes 
        std::vector<facet<3>> volumes = get_facets<3>(tet_a_t);
        ASSERT_EQ(volumes.size(), 1);
        ASSERT_EQ(volumes[0].t, tcode<3>{"000"});
        ASSERT_EQ(volumes[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2, 3}));
    }
    {
        // vertices
        std::vector<facet<0>> vertices = get_facets<0>(tet_b_t);
        ASSERT_EQ(vertices.size(), 4);
        ASSERT_EQ(vertices[0].t, tcode<0>{});
        ASSERT_EQ(vertices[0].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[0].vertex_indices, std::vector<std::size_t>{0});
        ASSERT_EQ(vertices[1].t, tcode<0>{});
        ASSERT_EQ(vertices[1].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[1].vertex_indices, std::vector<std::size_t>{1});
        ASSERT_EQ(vertices[2].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[2].vertex_indices, std::vector<std::size_t>{2});
        ASSERT_EQ(vertices[3].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[3].vertex_indices, std::vector<std::size_t>{3});

        // lines
        std::vector<facet<1>> lines = get_facets<1>(tet_b_t);
        ASSERT_EQ(lines.size(), 6);
        ASSERT_EQ(lines[0].t, tcode<1>{"1"});
        ASSERT_EQ(lines[0].vertex_indices, (std::vector<std::size_t>{0, 1}));
        ASSERT_EQ(lines[1].t, tcode<1>{"0"});
        ASSERT_EQ(lines[1].vertex_indices, (std::vector<std::size_t>{0, 2}));
        ASSERT_EQ(lines[2].t, tcode<1>{"0"});
        ASSERT_EQ(lines[2].vertex_indices, (std::vector<std::size_t>{1, 2}));
        ASSERT_EQ(lines[3].t, tcode<1>{"0"});
        ASSERT_EQ(lines[3].vertex_indices, (std::vector<std::size_t>{0, 3}));
        ASSERT_EQ(lines[4].t, tcode<1>{"0"});
        ASSERT_EQ(lines[4].vertex_indices, (std::vector<std::size_t>{1, 3}));
        ASSERT_EQ(lines[5].t, tcode<1>{"0"});
        ASSERT_EQ(lines[5].vertex_indices, (std::vector<std::size_t>{2, 3}));

        // surfaces
        std::vector<facet<2>> surfaces = get_facets<2>(tet_b_t);
        ASSERT_EQ(surfaces.size(), 4);
        ASSERT_EQ(surfaces[0].t, tcode<2>{"01"});
        ASSERT_EQ(surfaces[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2}));
        ASSERT_EQ(surfaces[1].t, tcode<2>{"01"});
        ASSERT_EQ(surfaces[1].vertex_indices, (std::vector<std::size_t>{0, 1, 3}));
        ASSERT_EQ(surfaces[2].t, tcode<2>{"00"});
        ASSERT_EQ(surfaces[2].vertex_indices, (std::vector<std::size_t>{0, 2, 3}));
        ASSERT_EQ(surfaces[3].t, tcode<2>{"00"});
        ASSERT_EQ(surfaces[3].vertex_indices, (std::vector<std::size_t>{1, 2, 3}));

        // volumes 
        std::vector<facet<3>> volumes = get_facets<3>(tet_b_t);
        ASSERT_EQ(volumes.size(), 1);
        ASSERT_EQ(volumes[0].t, tcode<3>{"001"});
        ASSERT_EQ(volumes[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2, 3}));
    }

    // ===========
    // = Pyramid =
    // ===========
    {
        // vertices
        std::vector<facet<0>> vertices = get_facets<0>(pyra_a_t);
        ASSERT_EQ(vertices.size(), 5);
        ASSERT_EQ(vertices[0].t, tcode<0>{});
        ASSERT_EQ(vertices[0].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[0].vertex_indices, std::vector<std::size_t>{0});
        ASSERT_EQ(vertices[1].t, tcode<0>{});
        ASSERT_EQ(vertices[1].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[1].vertex_indices, std::vector<std::size_t>{1});
        ASSERT_EQ(vertices[2].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[2].vertex_indices, std::vector<std::size_t>{2});
        ASSERT_EQ(vertices[3].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[3].vertex_indices, std::vector<std::size_t>{3});
        ASSERT_EQ(vertices[4].vertex_indices.size(), 1);
        ASSERT_EQ(vertices[4].vertex_indices, std::vector<std::size_t>{4});

        // lines
        std::vector<facet<1>> lines = get_facets<1>(pyra_a_t);
        ASSERT_EQ(lines.size(), 8);
        ASSERT_EQ(lines[0].t, tcode<1>{"0"});
        ASSERT_EQ(lines[0].vertex_indices, (std::vector<std::size_t>{0, 1}));
        ASSERT_EQ(lines[1].t, tcode<1>{"0"});
        ASSERT_EQ(lines[1].vertex_indices, (std::vector<std::size_t>{2, 3}));
        ASSERT_EQ(lines[2].t, tcode<1>{"1"});
        ASSERT_EQ(lines[2].vertex_indices, (std::vector<std::size_t>{0, 2}));
        ASSERT_EQ(lines[3].t, tcode<1>{"1"});
        ASSERT_EQ(lines[3].vertex_indices, (std::vector<std::size_t>{1, 3}));
        ASSERT_EQ(lines[4].t, tcode<1>{"0"});
        ASSERT_EQ(lines[4].vertex_indices, (std::vector<std::size_t>{0, 4}));
        ASSERT_EQ(lines[5].t, tcode<1>{"0"});
        ASSERT_EQ(lines[5].vertex_indices, (std::vector<std::size_t>{1, 4}));
        ASSERT_EQ(lines[6].t, tcode<1>{"0"});
        ASSERT_EQ(lines[6].vertex_indices, (std::vector<std::size_t>{2, 4}));
        ASSERT_EQ(lines[7].t, tcode<1>{"0"});
        ASSERT_EQ(lines[7].vertex_indices, (std::vector<std::size_t>{3, 4}));

        // surfaces
        std::vector<facet<2>> surfaces = get_facets<2>(pyra_a_t);
        ASSERT_EQ(surfaces.size(), 5);
        ASSERT_EQ(surfaces[0].t, tcode<2>{"10"});
        ASSERT_EQ(surfaces[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2, 3}));
        ASSERT_EQ(surfaces[1].t, tcode<2>{"00"});
        ASSERT_EQ(surfaces[1].vertex_indices, (std::vector<std::size_t>{0, 1, 4}));
        ASSERT_EQ(surfaces[2].t, tcode<2>{"00"});
        ASSERT_EQ(surfaces[2].vertex_indices, (std::vector<std::size_t>{2, 3, 4}));
        ASSERT_EQ(surfaces[3].t, tcode<2>{"01"});
        ASSERT_EQ(surfaces[3].vertex_indices, (std::vector<std::size_t>{0, 2, 4}));
        ASSERT_EQ(surfaces[4].t, tcode<2>{"01"});
        ASSERT_EQ(surfaces[4].vertex_indices, (std::vector<std::size_t>{1, 3, 4}));

        // volumes 
        std::vector<facet<3>> volumes = get_facets<3>(pyra_a_t);
        ASSERT_EQ(volumes.size(), 1);
        ASSERT_EQ(volumes[0].t, tcode<3>{"010"});
        ASSERT_EQ(volumes[0].vertex_indices, (std::vector<std::size_t>{0, 1, 2, 3, 4}));
    }
}

// Test that all "faces" (ndim - 1) have outwards unit normals
TEST(test_polytope, test_facet_normal){
    using namespace NUMTOOL::TENSOR::FIXED_SIZE;

    auto compute_centroid = []<geo_code auto t>(static_geo_code<t>){
        static constexpr std::size_t ndim = get_ndim(t);
        Tensor<double, ndim> centroid;
        centroid = 0.0;

        // contribution from each vertex
        for(vcode<ndim> v : gen_vert<t>()){
            for(int idim = 0; idim < ndim; ++idim){
                if(v[idim] != 0){
                    centroid[idim] += 1.0;
                } 
            }
        }

        // average
        for(int idim = 0; idim < ndim; ++idim){
            centroid[idim] /= n_vert(t);
        }

        return centroid;
    };

    auto compute_centroid_verts = []<geo_code auto t>(static_geo_code<t>,
            std::vector<std::size_t> verts){
        static constexpr std::size_t ndim = get_ndim(t);
        auto vert_list = gen_vert<t>();
        Tensor<double, ndim> centroid;
        centroid = 0.0;
        for(std::size_t ivert : verts){
            vcode<ndim> v = vert_list[ivert];
            for(int idim = 0; idim < ndim; ++idim){
                if(v[idim] != 0){
                    centroid[idim] += 1.0;
                } 
            }
        }
        // average
        for(int idim = 0; idim < ndim; ++idim){
            centroid[idim] /= verts.size();
        }

        return centroid;
    };

    auto check_normals = [=]<geo_code auto t>(static_geo_code<t> targ){
        static constexpr std::size_t ndim = get_ndim(t);

        Tensor<double, ndim> centroid = compute_centroid(targ);

        auto vert_list = gen_vert<t>();
        std::vector<facet<ndim - 1>> faces = get_facets<ndim - 1>(t);
        for(facet<ndim - 1> face : faces){
            Tensor<double, ndim, ndim - 1> axis_vectors;
            auto axis_nodes = orient_axis_vert(face.t);

            // get distance by subtract origin of orient axis 
            // and add each endpoint
            for(int j = 0; j < ndim - 1; ++j){
                std::size_t axis_v_idx = face.vertex_indices[axis_nodes[j+1]];
                std::size_t origin_v_idx = face.vertex_indices[axis_nodes[0]];
                auto axis_v = vert_list[axis_v_idx];
                auto origin_v = vert_list[origin_v_idx];
                for(int i = 0; i < ndim; ++i){
                    double axis_coord = axis_v[i] == 0 ? 0.0 : 1.0;
                    double origin_coord = origin_v[i] == 0 ? 0.0 : 1.0;
                    axis_vectors[i][j] = axis_coord - origin_coord;
                }
            }

            auto centroid_face = compute_centroid_verts(targ, face.vertex_indices);

            Tensor<double, ndim> normal{calc_ortho(axis_vectors)};

            
            SCOPED_TRACE("face vertex indices: "
                    + fmt::format("{}", face.vertex_indices) );
            SCOPED_TRACE("normal vector: "
                    + fmt::format("{}", normal) );
            ASSERT_TRUE(dot(normal, centroid_face - centroid) > 0);
        }
    };

    check_normals(static_geo_code<tri_a_t>{});

}

