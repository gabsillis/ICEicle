#include "Numtool/matrix/dense_matrix.hpp"
#include "Numtool/polydefs/LagrangePoly.hpp"
#include "iceicle/fe_function/nodal_fe_function.hpp"
#include <iceicle/transformations/HypercubeElementTransformation.hpp>
#include "gtest/gtest.h"
#include <random>
#include <limits>

using namespace ELEMENT::TRANSFORMATIONS;

TEST(test_hypercube_transform, test_ijk_poin){
  HypercubeElementTransformation<double, int, 1, 4> trans1{};
  ASSERT_EQ(
      "[ 0 ]\n"
      "[ 1 ]\n"
      "[ 2 ]\n"
      "[ 3 ]\n"
      "[ 4 ]\n",
      trans1.print_ijk_poin()
      );

  HypercubeElementTransformation<double, int, 3, 3> trans2{};
  ASSERT_EQ(
      "[ 0 0 0 ]\n"
      "[ 0 0 1 ]\n"
      "[ 0 0 2 ]\n"
      "[ 0 0 3 ]\n"
      "[ 0 1 0 ]\n"
      "[ 0 1 1 ]\n"
      "[ 0 1 2 ]\n"
      "[ 0 1 3 ]\n"
      "[ 0 2 0 ]\n"
      "[ 0 2 1 ]\n"
      "[ 0 2 2 ]\n"
      "[ 0 2 3 ]\n"
      "[ 0 3 0 ]\n"
      "[ 0 3 1 ]\n"
      "[ 0 3 2 ]\n"
      "[ 0 3 3 ]\n"
      "[ 1 0 0 ]\n"
      "[ 1 0 1 ]\n"
      "[ 1 0 2 ]\n"
      "[ 1 0 3 ]\n"
      "[ 1 1 0 ]\n"
      "[ 1 1 1 ]\n"
      "[ 1 1 2 ]\n"
      "[ 1 1 3 ]\n"
      "[ 1 2 0 ]\n"
      "[ 1 2 1 ]\n"
      "[ 1 2 2 ]\n"
      "[ 1 2 3 ]\n"
      "[ 1 3 0 ]\n"
      "[ 1 3 1 ]\n"
      "[ 1 3 2 ]\n"
      "[ 1 3 3 ]\n"
      "[ 2 0 0 ]\n"
      "[ 2 0 1 ]\n"
      "[ 2 0 2 ]\n"
      "[ 2 0 3 ]\n"
      "[ 2 1 0 ]\n"
      "[ 2 1 1 ]\n"
      "[ 2 1 2 ]\n"
      "[ 2 1 3 ]\n"
      "[ 2 2 0 ]\n"
      "[ 2 2 1 ]\n"
      "[ 2 2 2 ]\n"
      "[ 2 2 3 ]\n"
      "[ 2 3 0 ]\n"
      "[ 2 3 1 ]\n"
      "[ 2 3 2 ]\n"
      "[ 2 3 3 ]\n"
      "[ 3 0 0 ]\n"
      "[ 3 0 1 ]\n"
      "[ 3 0 2 ]\n"
      "[ 3 0 3 ]\n"
      "[ 3 1 0 ]\n"
      "[ 3 1 1 ]\n"
      "[ 3 1 2 ]\n"
      "[ 3 1 3 ]\n"
      "[ 3 2 0 ]\n"
      "[ 3 2 1 ]\n"
      "[ 3 2 2 ]\n"
      "[ 3 2 3 ]\n"
      "[ 3 3 0 ]\n"
      "[ 3 3 1 ]\n"
      "[ 3 3 2 ]\n"
      "[ 3 3 3 ]\n"
      , trans2.print_ijk_poin()
      );
}

TEST(test_hypercube_transform, test_fill_shp){
  static constexpr int Pn = 8;
  HypercubeElementTransformation<double, int, 4, Pn> trans1{};
  constexpr int nnode1 = trans1.n_nodes();
  std::array<double, nnode1> shp;
  MATH::GEOMETRY::Point<double, 4> xi = {0.3, 0.2, 0.1, 0.4};
  trans1.fill_shp(xi, shp.data());
  for(int inode = 0; inode < nnode1; ++inode){
    ASSERT_NEAR(
        (POLYNOMIAL::lagrange1d<double, Pn>(trans1.ijk_poin[inode][0], xi[0]) 
        * POLYNOMIAL::lagrange1d<double, Pn>(trans1.ijk_poin[inode][1], xi[1])  
        * POLYNOMIAL::lagrange1d<double, Pn>(trans1.ijk_poin[inode][2], xi[2])  
        * POLYNOMIAL::lagrange1d<double, Pn>(trans1.ijk_poin[inode][3], xi[3])),
      shp[inode], 1e-13
    );
  }
}

TEST( test_hypercube_transform, test_ref_coordinates ){
  HypercubeElementTransformation<double, int, 3, 2> trans1{};

  const MATH::GEOMETRY::Point<double, 3> *xi1 = trans1.reference_nodes();
  ASSERT_DOUBLE_EQ(-1.0, xi1[0][0]);
  ASSERT_DOUBLE_EQ(-1.0, xi1[0][1]);
  ASSERT_DOUBLE_EQ(-1.0, xi1[0][2]);

  ASSERT_DOUBLE_EQ(-1.0, xi1[1][0]);
  ASSERT_DOUBLE_EQ(-1.0, xi1[1][1]);
  ASSERT_DOUBLE_EQ( 0.0, xi1[1][2]);

  ASSERT_DOUBLE_EQ(-1.0, xi1[2][0]);
  ASSERT_DOUBLE_EQ(-1.0, xi1[2][1]);
  ASSERT_DOUBLE_EQ( 1.0, xi1[2][2]);

  ASSERT_DOUBLE_EQ(-1.0, xi1[3][0]);
  ASSERT_DOUBLE_EQ( 0.0, xi1[3][1]);
  ASSERT_DOUBLE_EQ(-1.0, xi1[3][2]);

  ASSERT_DOUBLE_EQ(-1.0, xi1[4][0]);
  ASSERT_DOUBLE_EQ( 0.0, xi1[4][1]);
  ASSERT_DOUBLE_EQ( 0.0, xi1[4][2]);

  ASSERT_DOUBLE_EQ(-1.0, xi1[5][0]);
  ASSERT_DOUBLE_EQ( 0.0, xi1[5][1]);
  ASSERT_DOUBLE_EQ( 1.0, xi1[5][2]);

  // TODO: check all
}

TEST( test_hypercube_transform, test_transform ){
  std::random_device rdev{};
  std::default_random_engine engine{rdev()};
  std::uniform_real_distribution<double> dist{-0.2, 0.2};
  std::uniform_real_distribution<double> domain_dist{-1.0, 1.0};

  auto rand_doub = [&]() -> double { return dist(engine); };

  HypercubeElementTransformation<double, int, 2, 1> trans_lin2d{};
  auto lagrange0 = [](double s){ return (1.0 - s) / 2.0; };
  auto lagrange1 = [](double s){ return (1.0 + s) / 2.0; };

  { // linear 2d transformation test
    int node_indices[trans_lin2d.n_nodes()];
    FE::NodalFEFunction<double, 2> node_coords(trans_lin2d.n_nodes());

    // randomly peturb
    for(int inode = 0; inode < trans_lin2d.n_nodes(); ++inode){
      node_indices[inode] = inode; // set the node index
      for(int idim = 0; idim < 2; ++idim){
        node_coords[inode][idim] = trans_lin2d.reference_nodes()[inode][idim] + rand_doub();
      }
    }
    for(int k = 0; k < 100; ++k){
      MATH::GEOMETRY::Point<double, 2> xi = {domain_dist(engine), domain_dist(engine)};
      MATH::GEOMETRY::Point<double, 2> x_act = {
        lagrange0(xi[0]) * lagrange0(xi[1]) * node_coords[0][0]
        + lagrange0(xi[0]) * lagrange1(xi[1]) * node_coords[1][0]
        + lagrange1(xi[0]) * lagrange0(xi[1]) * node_coords[2][0]
        + lagrange1(xi[0]) * lagrange1(xi[1]) * node_coords[3][0],

        lagrange0(xi[0]) * lagrange0(xi[1]) * node_coords[0][1]
        + lagrange0(xi[0]) * lagrange1(xi[1]) * node_coords[1][1]
        + lagrange1(xi[0]) * lagrange0(xi[1]) * node_coords[2][1]
        + lagrange1(xi[0]) * lagrange1(xi[1]) * node_coords[3][1]
      };

      MATH::GEOMETRY::Point<double, 2> x_trans;
      trans_lin2d.transform(node_coords, node_indices, xi, x_trans);
      ASSERT_DOUBLE_EQ(x_trans[0], x_act[0]);
      ASSERT_DOUBLE_EQ(x_trans[1], x_act[1]);
    }
  }

  HypercubeElementTransformation<double, int, 2, 3> trans1{};

  int node_indices[trans1.n_nodes()];
  FE::NodalFEFunction<double, 2> node_coords(trans1.n_nodes());

  // randomly peturb
  for(int inode = 0; inode < trans1.n_nodes(); ++inode){
    node_indices[inode] = inode; // set the node index
    for(int idim = 0; idim < 2; ++idim){
      node_coords[inode][idim] = trans1.reference_nodes()[inode][idim] + rand_doub();
    }
  }

  // Require: kronecker property
  for(int inode = 0; inode < trans1.n_nodes(); ++inode){
    const MATH::GEOMETRY::Point<double, 2> &xi = trans1.reference_nodes()[inode];
    MATH::GEOMETRY::Point<double, 2> x;
    trans1.transform(node_coords, node_indices, xi, x);
    for(int idim = 0; idim < 2; ++idim){
      ASSERT_NEAR(node_coords[inode][idim], x[idim], 1e-14); // tiny bit of roundoff error
    }
  }
}

TEST( test_hypercube_transform, test_fill_deriv ){
  static constexpr int ndim = 2;
  static constexpr int Pn = 1;
  HypercubeElementTransformation<double, int, ndim, Pn> trans{};

  auto lagrange0 = [](double s){ return (1.0 - s) / 2.0; };
  auto lagrange1 = [](double s){ return (1.0 + s) / 2.0; };

  auto dlagrange0 = [](double s){ return -0.5; };
  auto dlagrange1 = [](double s){ return 0.5; };

  MATH::MATRIX::DenseMatrixSetWidth<double, 2> dBidxj(4);
  MATH::MATRIX::DenseMatrixSetWidth<double, 2> dBidxj_2(4);
  MATH::GEOMETRY::Point<double, 2> xi = {0.3, -0.3};
  trans.fill_deriv(xi, dBidxj);
  trans.fill_deriv(xi, dBidxj_2);

  for(int inode = 0; inode < 4; ++inode) for(int idim = 0; idim < 2; ++idim) ASSERT_DOUBLE_EQ(dBidxj[inode][idim], dBidxj_2[inode][idim]);

  ASSERT_DOUBLE_EQ( dlagrange0(xi[0]) *  lagrange0(xi[1]), dBidxj[0][0]);
  ASSERT_DOUBLE_EQ(  lagrange0(xi[0]) * dlagrange0(xi[1]), dBidxj[0][1]);

  ASSERT_DOUBLE_EQ( dlagrange0(xi[0]) *  lagrange1(xi[1]), dBidxj[1][0]);
  ASSERT_DOUBLE_EQ(  lagrange0(xi[0]) * dlagrange1(xi[1]), dBidxj[1][1]);

  ASSERT_DOUBLE_EQ( dlagrange1(xi[0]) *  lagrange0(xi[1]), dBidxj[2][0]);
  ASSERT_DOUBLE_EQ(  lagrange1(xi[0]) * dlagrange0(xi[1]), dBidxj[2][1]);

  ASSERT_DOUBLE_EQ( dlagrange1(xi[0]) *  lagrange1(xi[1]), dBidxj[3][0]);
  ASSERT_DOUBLE_EQ(  lagrange1(xi[0]) * dlagrange1(xi[1]), dBidxj[3][1]);
}


TEST( test_hypercube_transform, test_jacobian ){
  std::random_device rdev{};
  std::default_random_engine engine{rdev()};
  std::uniform_real_distribution<double> dist{-0.2, 0.2};
  std::uniform_real_distribution<double> domain_dist{-1.0, 1.0};

  static double epsilon = 1e-8;//std::sqrt(std::numeric_limits<double>::epsilon());

  auto rand_doub = [&]() -> double { return dist(engine); };

  NUMTOOL::TMP::constexpr_for_range<1, 5>([&]<int ndim>(){
    NUMTOOL::TMP::constexpr_for_range<1, 9>([&]<int Pn>(){
      std::cout << "ndim: " << ndim << " | Pn: " << Pn << std::endl;
      HypercubeElementTransformation<double, int, ndim, Pn> trans1{};

      int node_indices[trans1.n_nodes()];
      FE::NodalFEFunction<double, ndim> node_coords(trans1.n_nodes());

      for(int k = 0; k < 50; ++k){ // repeat test 50 times
        // randomly peturb
        for(int inode = 0; inode < trans1.n_nodes(); ++inode){
          node_indices[inode] = inode; // set the node index
          for(int idim = 0; idim < ndim; ++idim){
            node_coords[inode][idim] = trans1.reference_nodes()[inode][idim] + rand_doub();
          }
        }

    //    std::cout << "pt1: [ " << node_coords[0][0] << ", " << node_coords[0][1] << " ]" << std::endl;
    //    std::cout << "pt2: [ " << node_coords[1][0] << ", " << node_coords[1][1] << " ]" << std::endl;
    //    std::cout << "pt3: [ " << node_coords[2][0] << ", " << node_coords[2][1] << " ]" << std::endl;
    //    std::cout << "pt4: [ " << node_coords[3][0] << ", " << node_coords[3][1] << " ]" << std::endl;
        // test 10 random points in the domain 
        for(int k2 = 0; k2 < 10; ++k2){
          MATH::GEOMETRY::Point<double, ndim> testpt;
          for(int idim = 0; idim < ndim; ++idim) testpt[idim] = domain_dist(engine);

          double Jtrans[ndim][ndim];
          double Jfd[ndim][ndim];
          trans1.Jacobian(node_coords, node_indices, testpt, Jtrans);

          MATH::GEOMETRY::Point<double, ndim> unpeturb_transform;
          trans1.transform(node_coords, node_indices, testpt, unpeturb_transform);

          for(int ixi = 0; ixi < ndim; ++ixi){
            double temp = testpt[ixi];
            testpt[ixi] += epsilon;
            MATH::GEOMETRY::Point<double, ndim> peturb_transform;
            trans1.transform(node_coords, node_indices, testpt, peturb_transform);
            for(int ix = 0; ix < ndim; ++ix){
              Jfd[ix][ixi] = (peturb_transform[ix] - unpeturb_transform[ix]) / epsilon;

            }
            testpt[ixi] = temp;
          }

    //      std::cout << Jtrans[0][0] << " " << Jtrans[0][1] << std::endl;
    //      std::cout << Jtrans[1][0] << " " << Jtrans[1][1] << std::endl;
    //      std::cout << std::endl;
    //
    //      std::cout << Jfd[0][0] << " " << Jfd[0][1] << std::endl;
    //      std::cout << Jfd[1][0] << " " << Jfd[1][1] << std::endl;
    //      std::cout << std::endl;
          double scaled_tol = 1e-6 * (std::pow(10, 0.4 * (Pn - 1)));
          for(int ix = 0; ix < ndim; ++ix) for(int ixi = 0; ixi < ndim; ++ixi) ASSERT_NEAR(Jtrans[ix][ixi], Jfd[ix][ixi], scaled_tol);
        }
      }
    });
  });
}

