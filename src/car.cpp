#include <Eigen/Dense>

#include <tassa_car.hpp>
#include <ilqg.hpp>

int main() {
  car::System sys;
  double dt = 0.03;
  const int max_iterations = 50;
  // const int T = 500;
  Vec4 x0 = (Vec4() << 1, 1, M_PI * .1, 0).finished();
  car::State state0(x0);
  car::State state_ref(Vec4::Zero());

  std::vector<car::U> list_u0(2000, .1 * Vec2::Ones());
  std::vector<car::State> list_ref(2000, state_ref);

  // Cost coefficients
  Mat4 Q = (Vec4() << .1, .1, 1, .3).finished().asDiagonal();
  Mat4 Qf = Q * 50;
  Mat2 R = (Vec2() << 1, .1).finished().asDiagonal();

  iLQG<car>::Params params(Q, R, Qf);

  Vec2 umin = (Vec2() << -.5, 2).finished();
  Vec2 umax = (Vec2() << .5, 2).finished();

  iLQG<car> ilqg(sys, dt);
  ilqg.init(state0, list_u0, list_ref, params, umin, umax, 10);
  ilqg.iterate(max_iterations, list_u0);
  return 0;
}
