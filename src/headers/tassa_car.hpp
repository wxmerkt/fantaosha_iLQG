#ifndef CAR
#define CAR

#include <list>
#include <iostream>
#include <cmath>
#include <string>
#include <functional>

#include <liegroup.hpp>

#include <mat.h>

#include <Eigen/Dense>

struct car {
  struct System;
  struct State;
  struct DState;

  static constexpr size_t M = 4;  // dimension of the configuration space
  static constexpr size_t N = 2;  // dimension of control inputs

  typedef Eigen::Matrix<double, M, 1> V;
  typedef Eigen::Matrix<double, N, 1> U;

  struct System {};

  //**************************************************************************************************************************
  // Definition of car state and dynamics
  // *************************************************************************************************************************

  struct State  // car state
      {
    Vec4 x;
    State(Vec4 x_ = Vec4::Zero());
    State update(const DState &dstate,
                 double h);  // compute next state by Euler method
    static Vec4 diff(State const &state1, State const &state2);
    static void save(const std::list<State> &, std::string);
  };

  struct DState  // time derivative of car state
      {
    Vec4 v;
    DState(const System &sys, const State &state,
           Vec2 u);  // compute body velocity
    DState(DState k1, DState k2, DState k3,
           DState k4);  // compute average body velocity for RK4
  };

  static void linearize(System const &sys, double const &dt, State const &state,
                        U const &u, Eigen::Matrix<double, M, M> &A,
                        Eigen::Matrix<double, M, N> &B) {
    double th = state.x(2);
    double v = state.x(3);

    double sth = sin(th);
    double cth = cos(th);

    A << 1, 0, -dt *sth *v, cth *dt, 0, 1, dt *cth *v, sth *dt, 0, 0, 1, 0, 0,
        0, 0, 1;

    // B = Eigen::Matrix<double, M, N>::Zero();
    // B(2, 0) = dt;
    // B(3, 1) = dt;
    B << 0, 0, 0, 0, dt, 0, 0, dt;
  }

  static double L(const Mat4 &M, const Mat2 &R, const Vec4 &dx, const Vec2 &u) {
    return (dx.transpose() * M * dx + u.transpose() * R * u)(0) * 0.5;
  }

  static Vec4 Lx(const Mat4 &M, const Vec4 &dx) { return M * dx; }

  static Vec2 Lu(const Mat2 &R, const Vec2 &u) { return R * u; }

  static Mat4 Lxx(const Mat4 &M, const Vec4 &dx) { return M; }

  static Mat2 Luu(const Mat2 &R, const Vec2 &u) { return R; }
};

car::State::State(Vec4 x_) : x(x_) {}

car::State car::State::update(const DState &dstate, double h) {
  State state_next;

  // TODO: investigate RK4
  state_next.x = this->x + dstate.v * h;

  return state_next;
}

Vec4 car::State::diff(State const &state1, State const &state2) {
  return state1.x - state2.x;
}
/*void car::State::save(const std::list<State> &states, std::string path) {
  MATFile *result;
  mxArray *thxy;
  void *p;
  mwSize *dims = new mwSize[2];

  result = matOpen(path.c_str(), "w");

  std::list<State>::const_iterator it_state;

  it_state = states.cbegin();
  dims[0] = 4;
  dims[1] = states.size();

  thxy = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

  for (p = mxGetPr(thxy); it_state != states.cend(); it_state++) {
    memcpy(p, it_state->x.data(), sizeof(double) * dims[0]);
    p += sizeof(double) * dims[0];
  }

  matPutVariable(result, "car", thxy);
  mxDestroyArray(thxy);
}*/

car::DState::DState(const System &sys, const State &state, Vec2 u) {
  // h-step dynamics are
  // x' = x + b cos theta
  // y' = y + b sin theta
  // theta' = theta + sin^-1 (sin(omega)f/d)
  // v' = v + ha

  Vec4 x = state.x;

  double d = 2.0;   // d = distance between back and front axles
  double h = 0.03;  // h = timestep (seconds)

  double w = u(0);  // w = front wheel angle
  double a = u(1);  // a = front wheel acceleration

  double o = x(2);  // o = car angle

  double v_ = x(3);   // v = front wheel velocity
  double f = h * v_;  // f = front wheel rolling distance
  double b =
      d + f * cos(w) -
      sqrt(d * d -
           (f * sin(w)) * (f * sin(w)));  // b = back wheel rolling distance
  double do_ = asin(sin(w) * f / d);      // do = change in car angle

  v << b *cos(o), b * sin(o), do_, h * a;
}

car::DState::DState(DState k1, DState k2, DState k3, DState k4) {
  std::cout << "EXCITING!" << std::endl;
  v = (k1.v + 2 * k2.v + 2 * k3.v + k4.v) / 6;
}
#endif
