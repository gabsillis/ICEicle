#include "iceicle/basis/lagrange.hpp"
#include "iceicle/element/reference_element.hpp"
#include "iceicle/fe_function/component_span.hpp"
#include "iceicle/fe_function/geo_layouts.hpp"
#include "iceicle/geometry/hypercube_element.hpp"
#include "iceicle/mesh/mesh.hpp"
#include "iceicle/tmp_utils.hpp"
#include <iceicle/quadrature/HypercubeGaussLegendre.hpp>
#include <iceicle/fe_function/fespan.hpp>
#include <iceicle/fe_function/dglayout.hpp>
#include <iceicle/fe_function/layout_right.hpp>
#include <iceicle/element/finite_element.hpp>

#include <gtest/gtest.h>

using namespace iceicle;
TEST(test_fespan, test_dglayout){

    using T = double;
    using IDX = int;
    static constexpr int ndim = 4;
    static constexpr int Pn = 3;
    using GeoElement = HypercubeElement<T, IDX, ndim, Pn>;
    using BasisType = HypercubeLagrangeBasis<T, IDX, ndim, Pn>;
    using QuadratureType = HypercubeGaussLegendre<T, IDX, ndim, Pn>;
    using FiniteElement = FiniteElement<T, IDX, ndim>;
    using EvaluationType = FEEvaluation<T, IDX, ndim>;

    std::vector<FiniteElement> elements;

    // create necesarry items to make finite elements 
    BasisType basis{};
    QuadratureType quadrule{};
    EvaluationType evals(basis, quadrule);

    // create geometric elements 
    GeoElement gel1{};
    GeoElement gel2{};
    int inode = 0;
    for(int i = 0; i < gel1.n_nodes(); ++i, ++inode){
        gel1.setNode(i, inode);
    }
    for(int i = 0; i < gel2.n_nodes(); ++i, ++inode){
        gel2.setNode(i, inode);
    }

    // create the finite elements 
    FiniteElement el1{&gel1, &basis, &quadrule, &evals, 0};
    FiniteElement el2{&gel2, &basis, &quadrule, &evals, 1};
    elements.push_back(el1);
    elements.push_back(el2);

    // get the offsets
    dg_dof_map offsets{elements};

    std::vector<T> data(offsets.calculate_size_requirement(2));
    std::iota(data.begin(), data.end(), 0.0);
    fe_layout_right<IDX, decltype(offsets), 2> layout(offsets);
    // alternate layout syntax
    using namespace tmp;
    fe_layout_right layout2{offsets, to_size<2>{}};

    fespan fespan1(data, layout);

    static constexpr int ndof_per_elem = ndim * (Pn + 1);
   
    int neq = 2;
    ASSERT_EQ(neq * 2 + 1.0, (fespan1[0, 2, 1]));

    std::size_t iel = 1, idof = 2, iv = 0;
    ASSERT_EQ(iel * std::pow(ndim, (Pn + 1)) * neq + idof * neq + iv, (fespan1[iel, idof, iv]));


    auto local_layout = fespan1.create_element_layout(1);
    std::vector<T> el_memory(local_layout.size());
    dofspan elspan1{el_memory.data(), local_layout};
    extract_elspan(1, fespan1, elspan1);

    iel = 1;
    idof = 8;
    iv = 1;
    ASSERT_DOUBLE_EQ(iel * std::pow(ndim, (Pn + 1)) * neq + idof * neq + iv, (elspan1[idof, iv]));

    scatter_elspan(1, 1.0, elspan1, 1.0, fespan1);
    ASSERT_DOUBLE_EQ(2.0 * (iel * std::pow(ndim, (Pn + 1)) * neq + idof * neq + iv), (fespan1[iel, idof, iv]));
    iel = 1;
    idof = 1;
    iv = 1;
    ASSERT_DOUBLE_EQ(2.0 *(iel * std::pow(ndim, (Pn + 1)) * neq + idof * neq + iv), (fespan1[iel, idof, iv]));
    iel = 0;
    idof = 8;
    iv = 0;
    // make sure original element data remains unchanged
    ASSERT_DOUBLE_EQ(iel * std::pow(ndim, (Pn + 1)) * neq + idof * neq + iv, (fespan1[iel, idof, iv]));

    //TODO: add more tests
}

TEST(test_dofspan, test_node_set_layout){
    using T = double;
    using IDX = int;
    static constexpr int ndim = 2;
    using namespace NUMTOOL::TENSOR::FIXED_SIZE;

    // set up an fespace with uniform mesh
    AbstractMesh<T, IDX, ndim> mesh{
        Tensor<T, 2>{{-1, -1}}, 
        Tensor<T, 2>{{1, 1}}, 
        Tensor<IDX, 2>{{4, 4}}
    };
    FESpace<T, IDX, ndim> fespace{&mesh, FESPACE_ENUMS::LAGRANGE, FESPACE_ENUMS::GAUSS_LEGENDRE, std::integral_constant<int, 1>{}};

    // Prerequisite: we assume a certain face->node connectivity
    // generated by the uniform mesh 
    
    ASSERT_EQ(mesh.faces[0]->nodes()[0], 1);
    ASSERT_EQ(mesh.faces[0]->nodes()[1], 6);

    ASSERT_EQ(mesh.faces[1]->nodes()[0], 2);
    ASSERT_EQ(mesh.faces[1]->nodes()[1], 7);

    ASSERT_EQ(mesh.faces[2]->nodes()[0], 3);
    ASSERT_EQ(mesh.faces[2]->nodes()[1], 8);

    ASSERT_EQ(mesh.faces[3]->nodes()[0], 6);
    ASSERT_EQ(mesh.faces[3]->nodes()[1], 11);

    ASSERT_EQ(mesh.faces[4]->nodes()[0], 7);
    ASSERT_EQ(mesh.faces[4]->nodes()[1], 12);

    ASSERT_EQ(mesh.faces[5]->nodes()[0], 8);
    ASSERT_EQ(mesh.faces[5]->nodes()[1], 13);

    ASSERT_EQ(mesh.faces[6]->nodes()[0], 11);
    ASSERT_EQ(mesh.faces[6]->nodes()[1], 16);

    ASSERT_EQ(mesh.faces[7]->nodes()[0], 12);
    ASSERT_EQ(mesh.faces[7]->nodes()[1], 17);

    ASSERT_EQ(mesh.faces[8]->nodes()[0], 13);
    ASSERT_EQ(mesh.faces[8]->nodes()[1], 18);

    ASSERT_EQ(mesh.faces[9]->nodes()[0], 16);
    ASSERT_EQ(mesh.faces[9]->nodes()[1], 21);

    ASSERT_EQ(mesh.faces[10]->nodes()[0], 17);
    ASSERT_EQ(mesh.faces[10]->nodes()[1], 22);

    ASSERT_EQ(mesh.faces[11]->nodes()[0], 18);
    ASSERT_EQ(mesh.faces[11]->nodes()[1], 23);

    ASSERT_EQ(mesh.faces[12]->nodes()[0], 6);
    ASSERT_EQ(mesh.faces[12]->nodes()[1], 5);

    ASSERT_EQ(mesh.faces[13]->nodes()[0], 7);
    ASSERT_EQ(mesh.faces[13]->nodes()[1], 6);

    ASSERT_EQ(mesh.faces[14]->nodes()[0], 8);
    ASSERT_EQ(mesh.faces[14]->nodes()[1], 7);

    ASSERT_EQ(mesh.faces[15]->nodes()[0], 9);
    ASSERT_EQ(mesh.faces[15]->nodes()[1], 8);

    ASSERT_EQ(mesh.faces[16]->nodes()[0], 11);
    ASSERT_EQ(mesh.faces[16]->nodes()[1], 10);

    ASSERT_EQ(mesh.faces[17]->nodes()[0], 12);
    ASSERT_EQ(mesh.faces[17]->nodes()[1], 11);

    ASSERT_EQ(mesh.faces[18]->nodes()[0], 13);
    ASSERT_EQ(mesh.faces[18]->nodes()[1], 12);

    ASSERT_EQ(mesh.faces[19]->nodes()[0], 14);
    ASSERT_EQ(mesh.faces[19]->nodes()[1], 13);

    ASSERT_EQ(mesh.faces[20]->nodes()[0], 16);
    ASSERT_EQ(mesh.faces[20]->nodes()[1], 15);

    ASSERT_EQ(mesh.faces[21]->nodes()[0], 17);
    ASSERT_EQ(mesh.faces[21]->nodes()[1], 16);

    ASSERT_EQ(mesh.faces[22]->nodes()[0], 18);
    ASSERT_EQ(mesh.faces[22]->nodes()[1], 17);

    ASSERT_EQ(mesh.faces[23]->nodes()[0], 19);
    ASSERT_EQ(mesh.faces[23]->nodes()[1], 18);

    ASSERT_EQ(mesh.faces[24]->nodes()[0], 5);
    ASSERT_EQ(mesh.faces[24]->nodes()[1], 0);

    ASSERT_EQ(mesh.faces[25]->nodes()[0], 4);
    ASSERT_EQ(mesh.faces[25]->nodes()[1], 9);

    ASSERT_EQ(mesh.faces[26]->nodes()[0], 10);
    ASSERT_EQ(mesh.faces[26]->nodes()[1], 5);

    ASSERT_EQ(mesh.faces[27]->nodes()[0], 9);
    ASSERT_EQ(mesh.faces[27]->nodes()[1], 14);

    ASSERT_EQ(mesh.faces[28]->nodes()[0], 15);
    ASSERT_EQ(mesh.faces[28]->nodes()[1], 10);

    ASSERT_EQ(mesh.faces[29]->nodes()[0], 14);
    ASSERT_EQ(mesh.faces[29]->nodes()[1], 19);

    ASSERT_EQ(mesh.faces[30]->nodes()[0], 20);
    ASSERT_EQ(mesh.faces[30]->nodes()[1], 15);

    ASSERT_EQ(mesh.faces[31]->nodes()[0], 19);
    ASSERT_EQ(mesh.faces[31]->nodes()[1], 24);

    ASSERT_EQ(mesh.faces[32]->nodes()[0], 0);
    ASSERT_EQ(mesh.faces[32]->nodes()[1], 1);

    ASSERT_EQ(mesh.faces[33]->nodes()[0], 21);
    ASSERT_EQ(mesh.faces[33]->nodes()[1], 20);

    ASSERT_EQ(mesh.faces[34]->nodes()[0], 1);
    ASSERT_EQ(mesh.faces[34]->nodes()[1], 2);

    ASSERT_EQ(mesh.faces[35]->nodes()[0], 22);
    ASSERT_EQ(mesh.faces[35]->nodes()[1], 21);

    ASSERT_EQ(mesh.faces[36]->nodes()[0], 2);
    ASSERT_EQ(mesh.faces[36]->nodes()[1], 3);

    ASSERT_EQ(mesh.faces[37]->nodes()[0], 23);
    ASSERT_EQ(mesh.faces[37]->nodes()[1], 22);

    ASSERT_EQ(mesh.faces[38]->nodes()[0], 3);
    ASSERT_EQ(mesh.faces[38]->nodes()[1], 4);

    ASSERT_EQ(mesh.faces[39]->nodes()[0], 24);
    ASSERT_EQ(mesh.faces[39]->nodes()[1], 23);


    // Try selecting some traces
    
    // =================
    //  Traces 4, 8, 21
    // =================

    nodeset_dof_map nodeset1{
        std::array<int, 3>{{4, 8, 21}},
        fespace
    };

    ASSERT_TRUE(std::ranges::contains(nodeset1.selected_traces, 4));
    ASSERT_TRUE(std::ranges::contains(nodeset1.selected_traces, 8));
    ASSERT_TRUE(std::ranges::contains(nodeset1.selected_traces, 21));
    ASSERT_EQ(nodeset1.selected_traces.size(), 3);

    // trace 4
    ASSERT_TRUE(std::ranges::contains(nodeset1.selected_nodes, 7));
    ASSERT_TRUE(std::ranges::contains(nodeset1.selected_nodes, 12));

    // trace 8
    ASSERT_TRUE(std::ranges::contains(nodeset1.selected_nodes, 13));
    ASSERT_TRUE(std::ranges::contains(nodeset1.selected_nodes, 18));

    // trace 21
    ASSERT_TRUE(std::ranges::contains(nodeset1.selected_nodes, 16));
    ASSERT_TRUE(std::ranges::contains(nodeset1.selected_nodes, 17));

    // inverse mapping 
    ASSERT_EQ(nodeset1.inv_selected_nodes[ 0], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[ 1], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[ 2], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[ 3], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[ 4], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[ 5], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[ 6], 6);
    ASSERT_TRUE(nodeset1.inv_selected_nodes[ 7] < 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[ 8], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[ 9], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[10], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[11], 6);
    ASSERT_TRUE(nodeset1.inv_selected_nodes[12] < 6);
    ASSERT_TRUE(nodeset1.inv_selected_nodes[13] < 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[14], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[15], 6);
    ASSERT_TRUE(nodeset1.inv_selected_nodes[16] < 6);
    ASSERT_TRUE(nodeset1.inv_selected_nodes[17] < 6);
    ASSERT_TRUE(nodeset1.inv_selected_nodes[18] < 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[19], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[20], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[21], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[22], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[23], 6);
    ASSERT_EQ(nodeset1.inv_selected_nodes[24], 6);

    // ========================
    //  Traces 5, 15, 17, 9, 0 
    // ========================
    nodeset_dof_map nodeset2{
        std::vector<int>{5, 15, 17, 9, 0},
        fespace
    };

    ASSERT_TRUE(std::ranges::contains(nodeset2.selected_traces, 5));
    ASSERT_TRUE(std::ranges::contains(nodeset2.selected_traces, 15));
    ASSERT_TRUE(std::ranges::contains(nodeset2.selected_traces, 17));
    ASSERT_TRUE(std::ranges::contains(nodeset2.selected_traces, 9));
    ASSERT_TRUE(std::ranges::contains(nodeset2.selected_traces, 0));

    // non boundary nodes 
    // 6, 8, 11, 12, 13, 16 
    ASSERT_EQ(nodeset2.selected_traces.size(), 5);
    ASSERT_EQ(nodeset2.selected_nodes.size(), 6);

    // inverse mapping 
    // NOTE: we don't explicitly state the ordering of selected_nodes 
    // but implementation results in in-order
    ASSERT_EQ(nodeset2.inv_selected_nodes[ 0], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[ 1], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[ 2], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[ 3], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[ 4], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[ 5], 6);
    ASSERT_TRUE(nodeset2.inv_selected_nodes[6] < 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[ 7], 6);
    ASSERT_TRUE(nodeset2.inv_selected_nodes[8] < 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[ 9], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[10], 6);
    ASSERT_TRUE(nodeset2.inv_selected_nodes[11] < 6);
    ASSERT_TRUE(nodeset2.inv_selected_nodes[12] < 6);
    ASSERT_TRUE(nodeset2.inv_selected_nodes[13] < 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[14], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[15], 6);
    ASSERT_TRUE(nodeset2.inv_selected_nodes[16] < 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[17], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[18], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[19], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[20], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[21], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[22], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[23], 6);
    ASSERT_EQ(nodeset2.inv_selected_nodes[24], 6);
}


TEST(test_dofspan, test_geo_dof_map){
    using T = double;
    using IDX = int;
    static constexpr int ndim = 2;
    using namespace NUMTOOL::TENSOR::FIXED_SIZE;

    // set up an fespace with uniform mesh
    AbstractMesh<T, IDX, ndim> mesh{
        Tensor<T, 2>{{-1, -1}}, 
        Tensor<T, 2>{{1, 1}}, 
        Tensor<IDX, 2>{{4, 4}}
    };
    FESpace<T, IDX, ndim> fespace{&mesh, FESPACE_ENUMS::LAGRANGE, FESPACE_ENUMS::GAUSS_LEGENDRE, std::integral_constant<int, 1>{}};

    fixed_component_constraint<T, ndim> left_wall_constraint{-1.0, 0};
    geo_dof_map geo_map{std::vector<int>{5, 15, 17, 9, 0}, fespace};

    ASSERT_TRUE(std::ranges::contains(geo_map.selected_traces, 5));
    ASSERT_TRUE(std::ranges::contains(geo_map.selected_traces, 15));
    ASSERT_TRUE(std::ranges::contains(geo_map.selected_traces, 17));
    ASSERT_TRUE(std::ranges::contains(geo_map.selected_traces, 9));
    ASSERT_TRUE(std::ranges::contains(geo_map.selected_traces, 0));

    // non boundary nodes 
    // 6, 8, 11, 12, 13, 16 
    ASSERT_EQ(geo_map.selected_traces.size(), 5);
    ASSERT_EQ(geo_map.selected_nodes.size(), 9);

    // inverse mapping 
    // NOTE: we don't explicitly state the ordering of selected_nodes 
    // but implementation results in in-order
    ASSERT_EQ(  geo_map.inv_selected_nodes[ 0], 9);
    ASSERT_TRUE(geo_map.inv_selected_nodes[ 1] < 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[ 2], 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[ 3], 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[ 4], 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[ 5], 9);
    ASSERT_TRUE(geo_map.inv_selected_nodes[ 6] < 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[ 7], 9);
    ASSERT_TRUE(geo_map.inv_selected_nodes[ 8] < 9);
    ASSERT_TRUE(geo_map.inv_selected_nodes[ 9] < 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[10], 9);
    ASSERT_TRUE(geo_map.inv_selected_nodes[11] < 9);
    ASSERT_TRUE(geo_map.inv_selected_nodes[12] < 9);
    ASSERT_TRUE(geo_map.inv_selected_nodes[13] < 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[14], 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[15], 9);
    ASSERT_TRUE(geo_map.inv_selected_nodes[16] < 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[17], 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[18], 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[19], 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[20], 9);
    ASSERT_TRUE(geo_map.inv_selected_nodes[21] < 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[22], 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[23], 9);
    ASSERT_EQ(  geo_map.inv_selected_nodes[24], 9);

    geo_data_layout layout{geo_map};
    std::vector<T> storage(layout.size());
    component_span data{storage, layout};
    extract_geospan(mesh, data);

    // set node 12 (0.0, 0.0) to something different
    int ldof12 = geo_map.inv_selected_nodes[12];
    ASSERT_DOUBLE_EQ((data[ldof12, 0]), 0.0);
    ASSERT_DOUBLE_EQ((data[ldof12, 1]), 0.0);
    data[ldof12, 1] = 0.05;

    update_mesh(data, mesh);
    ASSERT_DOUBLE_EQ(mesh.nodes[12][1], 0.05);

}
