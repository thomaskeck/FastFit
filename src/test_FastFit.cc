/**
 * Thomas Keck 2016
 */

#include <gtest/gtest.h>
#include <FastFit.h>

/** Test fixture. */
class FastFitTest : public ::testing::Test {
protected:
  virtual void SetUp()
  {
  }

  virtual void TearDown()
  {
  }
};


TEST_F(FastFitTest, HelixWithoutMagneticField) {

    // See HelixParametrisation.ipynb for an independent calculation of the parameters using sympy
    Eigen::Matrix<double, 3, 1> x;
    x << 3.0, 2.0, 1.0;  
    Eigen::Matrix<double, 3, 1> p;
    p << 1.0, 2.0, 3.0;  

    // We move the position each iteration along the momentum direction
    // to test if the parametrisation is stable
    for(unsigned int i = 0; i < 10; ++i) {

        Helix helix(0.0, x, p);
        auto m = helix.GetParametrisation();
        EXPECT_NEAR(m(0),  0.44721, 0.001);
        EXPECT_NEAR(m(1),  1.34164, 0.001);
        EXPECT_NEAR(m(2),  1.10714, 0.001);
        EXPECT_NEAR(m(3), -1.10557, 0.001);
        EXPECT_NEAR(m(4), -3.20000, 0.001);
        
        x += p;
    }

    // Reset position because some elements of the Jacobian are not invariant!
    x << 3.0, 2.0, 1.0;  
    Helix helix(0.0, x, p);
    auto A = helix.GetVertexJacobian();
    std::vector<std::vector<float>> expected_A = {{0, 0, 0},
                                                  {0, 0, 0},
                                                  {0, 0, 0},
                                                  {-0.552786404500042, 0.276393202250021, 0},
                                                  {-0.600000000000000, -1.20000000000000, 1}};
    for(unsigned int i = 0; i < 5; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
            EXPECT_NEAR(A(i,j), expected_A[i][j], 0.001);
        }
    }
    
    auto B = helix.GetMomentumJacobian();
    std::vector<std::vector<float>> expected_B = {{-0.0894427190999916, -0.178885438199983, 0},
                                                  {-0.268328157299975, -0.536656314599950, 0.447213595499958},
                                                  {-0.400000000000000, 0.200000000000000, 0},
                                                  {0.926687370800101, -0.0813776741499454, 0},
                                                  {-0.120000000000000, 2.16000000000000, -1.4}};
    for(unsigned int i = 0; i < 5; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
            EXPECT_NEAR(B(i,j), expected_B[i][j], 0.001);
        }
    }

}

TEST_F(FastFitTest, HelixWithMagneticField) {

    // See HelixParametrisation.ipynb for an independent calculation of the parameters using sympy
    Eigen::Matrix<double, 3, 1> x;
    x << 3.0, 2.0, 1.0;  
    Eigen::Matrix<double, 3, 1> p;
    p << 1.0, 2.0, 3.0;  

    const double pt = std::sqrt(p(0)*p(0) + p(1)*p(1));

    // We move the position and momentum each iteration along the helix
    // to test if the parametrisation is stable
    for(unsigned int i = 0; i < 10; ++i) {

        const double phi = 0.1 * static_cast<float>(i);
        p(0) = -std::sin(phi) * pt;
        p(1) =  std::cos(phi) * pt;
        x(0) = 3.0 + std::cos(phi) * pt;
        x(1) = 2.0 + std::sin(phi) * pt;
        x(2) = 1.0 + p(2) * phi;

        Helix helix(1.0, x, p);
        auto m = helix.GetParametrisation();
        EXPECT_NEAR(m(0),  0.44721, 0.001);
        EXPECT_NEAR(m(1),  1.34164, 0.001);
        EXPECT_NEAR(m(2), -0.98279, 0.001);
        EXPECT_NEAR(m(3),  0.92963, 0.001);

        // There is a jump in the z0 component at one point
        // Jumps are fine, because the parametrisation is not unique if the helix describes more than one circle in the considered region,
        // in the detector this is usually never the case (it could happen for curlers though)
        if( m(4) < 0 ) {
            EXPECT_NEAR(m(4), -6.66077, 0.001);
        } else {
            EXPECT_NEAR(m(4), 12.18878, 0.001);
        }
        
        x += p;
    }
    
    // Reset position and momentum because some elements of the Jacobian are not invariant!
    x << 3.0, 2.0, 1.0;  
    p << 1.0, 2.0, 3.0;
    const double phi = 0; 
    p(0) = -std::sin(phi) * pt;
    p(1) =  std::cos(phi) * pt;
    x(0) = 3.0 + std::cos(phi) * pt;
    x(1) = 2.0 + std::sin(phi) * pt;
    x(2) = 1.0 + p(2) * phi;
    Helix helix(1.0, x, p);
    auto A = helix.GetVertexJacobian();
    std::vector<std::vector<float>> expected_A = {{0, 0, 0},
                                                  {0, 0, 0},
                                                  {-0.153846153846154, 0.230769230769231, 0},
                                                  {0.607340407554780, 0.404893605036520, 0},
                                                  {-0.461538461538462, 0.692307692307692, 1}};
    for(unsigned int i = 0; i < 5; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
            EXPECT_NEAR(A(i,j), expected_A[i][j], 0.001);
        }
    }
    
    auto B = helix.GetMomentumJacobian();
    std::vector<std::vector<float>> expected_B = {{0, -0.200000000000000, 0},
                                                  {0, -0.600000000000000, 0.447213595499958},
                                                  {0.230769230769231, 0.153846153846154, 0},
                                                  {0.404893605036520, -1.61013250431113, 0},
                                                  {2.03394847880757, 0.461538461538462, -2.55359005004}};
    for(unsigned int i = 0; i < 5; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
            EXPECT_NEAR(B(i,j), expected_B[i][j], 0.001);
        }
    }
}

std::vector<std::vector<double>> createDiagonalErrorMatrix() {
  
  std::vector<std::vector<double>> error;
  for(unsigned int i = 0; i < 7; ++i) {
    std::vector<double> temp(7);
    for(unsigned int j = 0; j < 7; ++j)
      temp[j] = (i == j) ? 0.01 : 0.0;
    error.push_back(temp);
  }

  return error;

}

TEST_F(FastFitTest, TestNeutralParticles)
{
  return;
  FastFit fitter(2);

  fitter.SetDaughter(0, 0, std::vector<double>{-1.0, 0.0, 1.0}, std::vector<double>{ 1.0, 0.0, 0.0}, createDiagonalErrorMatrix());
  fitter.SetDaughter(1, 0, std::vector<double>{ 1.0, 0.0, 1.0}, std::vector<double>{-1.0, 0.0, 0.0}, createDiagonalErrorMatrix());
  fitter.fit(3, 1.5);
  
  // Neutral particles like photons just move along a straight line
  // So The solution is just the intersection point of the two lines

  EXPECT_NEAR(fitter.GetVertex(0), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(1), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(2), 1.0, 0.01);
  
  // Momenta should be the same

  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 0), -1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 1),  0.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 2),  1.0, 0.01);
  
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 0),  1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 1),  0.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 2),  1.0, 0.01);
  
  // Neutral particles like photons just move along a straight line
  // So The solution is just the intersection point of the two lines
  
  fitter.SetDaughter(0, 0, std::vector<double>{0.0, 1.0, -1.0}, std::vector<double>{0.0, 0.0,  1.0}, createDiagonalErrorMatrix());
  fitter.SetDaughter(1, 0, std::vector<double>{0.0, 1.0,  1.0}, std::vector<double>{0.0, 0.0, -1.0}, createDiagonalErrorMatrix());
  fitter.fit(3, 1.5);

  EXPECT_NEAR(fitter.GetVertex(0), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(1), 1.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(2), 0.0, 0.01);
  
  // Momenta should be the same

  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 0),  0.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 1),  1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 2), -1.0, 0.01);
  
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 0),  0.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 1),  1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 2),  1.0, 0.01);
  
  // Neutral particles like photons just move along a straight line
  // So The solution is just the intersection point of the two lines
  
  fitter.SetDaughter(0, 0, std::vector<double>{1.0, -1.0, 1.0}, std::vector<double>{0.0,  1.0, 0.0}, createDiagonalErrorMatrix());
  fitter.SetDaughter(1, 0, std::vector<double>{1.0,  1.0, 1.0}, std::vector<double>{0.0, -1.0, 0.0}, createDiagonalErrorMatrix());
  fitter.fit(3, 1.5);

  EXPECT_NEAR(fitter.GetVertex(0), 1.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(1), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(2), 1.0, 0.01);
  
  // Momenta should be the same

  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 0),  1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 1), -1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 2),  1.0, 0.01);
  
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 0),  1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 1),  1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 2),  1.0, 0.01);

  // Now what happens if the momenta do not match?
  // In the following case y momentum should to -0.1 and 0.1
  // In this case the momenta of the daughter particles have to be corrected
  
  fitter.SetDaughter(0, 0, std::vector<double>{-1.0, 0.0, 1.0}, std::vector<double>{ 1.0,  0.1, 0.0}, createDiagonalErrorMatrix());
  fitter.SetDaughter(1, 0, std::vector<double>{ 1.0, 0.0, 1.0}, std::vector<double>{-1.0, -0.1, 0.0}, createDiagonalErrorMatrix());
  fitter.fit(3, 1.5);

  EXPECT_NEAR(fitter.GetVertex(0), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(1), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(2), 1.0, 0.01);
  
  // Momenta should adapted so that they match

  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 0), -1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 1), -0.1, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(0, 2),  1.0, 0.01);
  
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 0),  1.0, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 1),  0.1, 0.01);
  EXPECT_NEAR(fitter.GetDaughterMomentum(1, 2),  1.0, 0.01);

}

TEST_F(FastFitTest, TestChargedParticles)
{
  return;
  FastFit fitter(2);

  for (unsigned int i = 1; i < 100; ++i) {
    fitter.SetDaughter(0,  1, std::vector<double>{0.0,  0.1 * i, 0.01}, std::vector<double>{0.0,  10.0, 10.0}, createDiagonalErrorMatrix());
    fitter.SetDaughter(1, -1, std::vector<double>{0.0, -0.1 * i, 0.01}, std::vector<double>{0.0, -10.0, 10.0}, createDiagonalErrorMatrix());
    fitter.fit(3, 1.5);

    //std::cout << "Momentum in y direction " << i * 0.1 << std::endl;
    //std::cout << fitter.GetVertex(0) << " " << fitter.GetVertex(1) << " " << fitter.GetVertex(2) << " " << std::endl;
  }

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
