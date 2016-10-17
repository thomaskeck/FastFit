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


TEST_F(FastFitTest, TestSimpleFit)
{

  FastFit fitter(2);

  std::vector<std::vector<double>> error;
  for(unsigned int i = 0; i < 7; ++i) {
    std::vector<double> temp(7);
    for(unsigned int j = 0; j < 7; ++j)
      temp[j] = (i == j) ? 1.0 : 0.0;
    error.push_back(temp);
  }

  fitter.SetDaughter(0, 0, std::vector<double>{-1.0, 0.0, 1.0}, std::vector<double>{ 1.0, 0.0, 0.0}, error);
  fitter.SetDaughter(1, 0, std::vector<double>{ 1.0, 0.0, 1.0}, std::vector<double>{-1.0, 0.0, 0.0}, error);
  fitter.fit(3, 1.5);

  EXPECT_NEAR(fitter.GetVertex(0), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(1), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(2), 1.0, 0.01);

  fitter.SetDaughter(0, 0, std::vector<double>{0.0, 1.0, -1.0}, std::vector<double>{0.0, 0.0,  1.0}, error);
  fitter.SetDaughter(1, 0, std::vector<double>{0.0, 0.0,  1.0}, std::vector<double>{0.0, 0.0, -1.0}, error);
  fitter.fit(3, 1.5);

  EXPECT_NEAR(fitter.GetVertex(0), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(1), 1.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(2), 0.0, 0.01);
  
  fitter.SetDaughter(0, 0, std::vector<double>{1.0, -1.0, 1.0}, std::vector<double>{0.0,  1.0, 0.0}, error);
  fitter.SetDaughter(1, 0, std::vector<double>{1.0,  1.0, 1.0}, std::vector<double>{0.0, -1.0, 0.0}, error);
  fitter.fit(3, 1.5);

  EXPECT_NEAR(fitter.GetVertex(0), 1.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(1), 0.0, 0.01);
  EXPECT_NEAR(fitter.GetVertex(2), 1.0, 0.01);


  for (unsigned int i = 1; i < 100; ++i) {
    fitter.SetDaughter(0,  1, std::vector<double>{0.0,  0.1 * i, 0.01}, std::vector<double>{0.0,  10.0, 10.0}, error);
    fitter.SetDaughter(1, -1, std::vector<double>{0.0, -0.1 * i, 0.01}, std::vector<double>{0.0, -10.0, 10.0}, error);
    fitter.fit(3, 1.5);

    std::cout << "Momentum in y direction " << i * 0.1 << std::endl;
    std::cout << fitter.GetVertex(0) << " " << fitter.GetVertex(1) << " " << fitter.GetVertex(2) << " " << std::endl;
  }

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
