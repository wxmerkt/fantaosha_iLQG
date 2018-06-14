#include <list>
#include <iostream>
#include <cmath>
#include <string>
#include <functional>

#include <mat.h>

#include <Eigen/Dense>

struct urdf_robot {
  struct System;
  struct State;
  struct DState;

  static constexpr size_t M = 5;  // dimension of the configuration space
  static constexpr size_t N = 2;  // dimension of control inputs

  typedef Eigen::Matrix<double, M, 1> V;
  typedef Eigen::Matrix<double, N, 1> U;

  struct System {};

  /**
   * @brief      URDF Robot State
   */
  struct State {
    Vec5 x;
    State(Vec5 x_ = Vec5::Zero());
    State update(const DState &dstate,
                 double h);  // compute next state by Euler method
    static Vec5 diff(State const &state1, State const &state2);
    static void save(const std::list<State> &, std::string);
  };

  /**
   * @brief      Time derivative of URDF Robot State
   */
  struct DState {
    Vec5 v;
    DState(const System &sys, const State &state,
           Vec2 u);  // compute body velocity
  };

  // TODO: LINEARIZE!!
  static void linearize(System const &sys, double const &dt, State const &state,
                        U const &u, Eigen::Matrix<double, 5, 5> &A,
                        Eigen::Matrix<double, 5, 2> &B) {
    double th = state.x(2);
    double v = state.x(3);

    double sth = sin(th);
    double cth = cos(th);

    A << 1, 0, -dt *sth *v, cth *dt, 0, 0, 1, dt *cth *v, sth *dt, 0, 0, 0, 1,
        0, dt, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1;

    B << 0, 0, 0, 0, 0, 0, dt, 0, 0, dt;
  }

  static double L(const Mat5 &M, const Mat2 &R, const Vec5 &dx, const Vec2 &u) {
    return (dx.transpose() * M * dx + u.transpose() * R * u)(0) * 0.5;
  }

  static Vec5 Lx(const Mat5 &M, const Vec5 &dx) { return M * dx; }

  static Vec2 Lu(const Mat2 &R, const Vec2 &u) { return R * u; }

  static Mat5 Lxx(const Mat5 &M, const Vec5 &dx) { return M; }

  static Mat2 Luu(const Mat2 &R, const Vec2 &u) { return R; }
};

urdf_robot::State::State(Vec5 x_) : x(x_) {}

urdf_robot::State urdf_robot::State::update(const DState &dstate, double h) {
  State state_next;

  state_next.x = this->x + dstate.v * h;

  return state_next;
}

Vec5 urdf_robot::State::diff(State const &state1, State const &state2) {
  return state1.x - state2.x;
}

/**
 * @brief      Saves robot state to Matlab file
 *
 * @param[in]  states  List of states to save (i.e. state trajectory)
 * @param[in]  path    Path of the Matlab output file
 */
void urdf_robot::State::save(const std::list<State> &states, std::string path) {
  MATFile *result;

  mxArray *thxy;

  void *p;

  mwSize *dims = new mwSize[2];

  result = matOpen(path.c_str(), "w");

  std::list<State>::const_iterator it_state;

  it_state = states.cbegin();
  dims[0] = 5;
  dims[1] = states.size();

  thxy = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

  for (p = mxGetPr(thxy); it_state != states.cend(); it_state++) {
    memcpy(p, it_state->x.data(), sizeof(double) * dims[0]);
    p += sizeof(double) * dims[0];
  }

  matPutVariable(result, "car", thxy);
  mxDestroyArray(thxy);
}

/**
 * @brief      Time derivative of Robot State
 *
 * @param[in]  sys    The system
 * @param[in]  state  The state
 * @param[in]  u      The applied controls
 */
urdf_robot::DState::DState(const System &sys, const State &state, Vec2 u) {
  double th = state.x(2);
  v << cos(th) * state.x(3), sin(th) * state.x(3), state.x(4), u;
}
