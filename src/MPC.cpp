#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;
using Eigen::VectorXd;

/**
 * TODO: Set the timestep length and duration
 */
size_t N = 25;
double dt = 0.08;

double ref_v = 40;
// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
//   simulator around in a circle with a constant steering angle and velocity on
//   a flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
//   presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

int x_start = 0;
int y_start = x_start + N;
int phi_start = y_start + N;
int v_start = phi_start + N;
int cte_y_start = v_start + N;
int cte_phi_start = cte_y_start + N;
int delta_start = cte_phi_start + N;
int a_start = delta_start + N - 1;

double v_ref = 40;
class FG_eval
{
public:
  // Fitted polynomial coefficients
  VectorXd coeffs;
  FG_eval(VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector &fg, const ADvector &vars)
  {
    /**
     * TODO: implement MPC
     * `fg` is a vector of the cost constraints, `vars` is a vector of variable 
     *   values (state & actuators)
     * NOTE: You'll probably go back and forth between this function and
     *   the Solver function below.
     */
    //Generell cost
    fg[0] = 0;
    for (int t = 0; t < N; ++t)
    { //Not being in lane centre
      fg[0] += CppAD::pow(vars[cte_y_start + t], 2);
      //Not being oriented to lane
      fg[0] += CppAD::pow(vars[cte_phi_start + t], 2);
      //Not being at reference speed
      fg[0] += CppAD::pow(vars[v_start + t] - v_ref, 2);
    }

    //Cost for using actuators (minimize energy)
    for (int t = 0; t < N - 1; ++t)
    { //For steering
      fg[0] += CppAD::pow(vars[delta_start + t], 2);
      //For accelerating/braking
      fg[0] += CppAD::pow(vars[a_start + t], 2);
    }

    //Cost for to high differences between actuations (smoother controlls)
    for (int t = 0; t < N - 2; ++t)
    { //For steering
      fg[0] += CppAD::pow(vars[delta_start + 1] - vars[delta_start + t], 2);
      //For accelerating/braking
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    //Init constraints
    fg[x_start + 1] = vars[x_start];
    fg[y_start + 1] = vars[y_start];
    fg[phi_start + 1] = vars[phi_start];
    fg[v_start + 1] = vars[v_start];
    fg[cte_y_start + 1] = vars[cte_y_start];
    fg[cte_phi_start + 1] = vars[cte_phi_start];

    //Rest of the constraints
    for (size_t t = 1; t < N; t++)
    {
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> phi1 = vars[phi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte_y1 = vars[cte_y_start + t];
      AD<double> cte_phi1 = vars[cte_phi_start + t];

      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> phi0 = vars[phi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte_y0 = vars[cte_y_start + t - 1];
      AD<double> cte_phi0 = vars[cte_phi_start + t - 1];

      AD<double> a0 = vars[a_start + t - 1];
      AD<double> delta0 = vars[delta_start + t - 1];

      fg[x_start + t + 1] = x1 - (x0 * CppAD::cos(phi0) * dt);
      fg[y_start + t + 1] = y1 - (y0 * CppAD::sin(phi0) * dt);
      fg[phi_start + t + 1] = phi1 - (phi0 + (v0 / Lf) * delta0 * dt);
      fg[x_start + t + 1] = v1 - (v0 + a0 * dt);

      AD<double> f0 = coeffs[1] * x0 + coeffs[0];
      fg[cte_y_start + t + 1] = cte_y1 - (f0 - y0 + v0 * CppAD::sin(cte_phi0));

      AD<double> dphi0 = CppAD::atan(coeffs[1]);
      fg[cte_phi_start + t + 1] = cte_phi1 - (phi0 - dphi0 + (v0 / Lf) * delta0 * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

std::vector<double> MPC::Solve(const VectorXd &state, const VectorXd &coeffs)
{
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  /**
   * TODO: Set the number of model variables (includes both states and inputs).
   * For example: If the state is a 4 element vector, the actuators is a 2
   *   element vector and there are 10 timesteps. The number of variables is:
   *   4 * 10 + 2 * 9
   */
  size_t n_vars = 4 * N + 2 * (N - 1);
  /**
   * TODO: Set the number of constraints
   */
  size_t n_constraints = 6 * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; ++i)
  {
    vars[i] = 0;
  }
  // Set the initial variable values
  vars[x_start] = state[0];
  vars[y_start] = state[1];
  vars[phi_start] = state[2];
  vars[v_start] = state[3];
  vars[cte_y_start] = state[4];
  vars[cte_phi_start] = state[5];

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  /**
   * TODO: Set lower and upper limits for variables.
   */

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; ++i)
  {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; ++i)
  {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; ++i)
  {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; ++i)
  {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  // NOTE: You don't have to worry about these options
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  //   of sparse routines, this makes the computation MUCH FASTER. If you can
  //   uncomment 1 of these and see if it makes a difference or not but if you
  //   uncomment both the computation time should go up in orders of magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  return {solution.x[x_start + 1], solution.x[y_start + 1],
          solution.x[phi_start + 1], solution.x[v_start + 1],
          solution.x[cte_y_start + 1], solution.x[cte_phi_start + 1],
          solution.x[delta_start], solution.x[a_start]};
  return {};
}