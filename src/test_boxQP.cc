#include <iostream>
#include <boxQP.hpp>

static constexpr int n = 100;
typedef Eigen::Matrix<double, n, 1> Vec;
typedef Eigen::Matrix<double, n, n> Mat;

int main() {
  // Uses uniform random numbers instead of Gaussian
  Vec g;
  g.setRandom();
  Mat H;
  H.setRandom();
  H = 5 * H;
  H = H * H.inverse();
  Vec lower = -Vec::Ones();
  Vec upper;
  upper.setOnes();

  Vec x0;
  x0.setRandom();

  Vec x;
  x.setZero();
  Eigen::MatrixXd Hfree;
  Eigen::Array<int, n, 1> free;
  int result = -1;

  boxQP<n>(H, g, lower, upper, x0, x, Hfree, free, result);
  std::cout << "Result: " << result << std::endl;
  return 0;
}