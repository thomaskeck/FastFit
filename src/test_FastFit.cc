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
