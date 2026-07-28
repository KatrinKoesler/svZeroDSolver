// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "ValveSpline.h"

void ValveSpline::setup_dofs(DOFHandler& dofhandler) {
  // set_up_dofs args: dofhandler (passed in), num equations, list of internal
  // variable names (strings) 3 eqns, one for Pressure, one for Flow, one for
  // the valve status output
  Block::setup_dofs_(dofhandler, 3, {"valve_status"});
}

// update_constant updates matrices E and F from E(y,t)*y_dot + F(y,t)*y +
// c(y,t) = 0 with terms that DO NOT DEPEND ON THE SOLUTION
void ValveSpline::update_constant(SparseSystem& system,
                                  std::vector<double>& parameters) {
  // Set element contributions
  // coeffRef args are the indices (i,j) of the matrix
  // global_eqn_ids: number of rows in the matrix, set in setup_dofs
  // global_var_ids: number of columns, organized as pressure and flow of all
  // inlets and then all outlets of the block
  double Rmin = parameters[global_param_ids[ParamId::RMIN]];
  double Rmax = parameters[global_param_ids[ParamId::RMAX]];
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[4]) = 1.0;
}

// update_solution updates matrices E and F from E(y,t)*y_dot + F(y,t)*y +
// c(y,t) = 0 with terms that DO DEPEND ON THE SOLUTION (will change with each
// time step)
void ValveSpline::update_solution(
    SparseSystem& system, std::vector<double>& parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) {
  // Get states
  double p_in = y[global_var_ids[0]];
  double p_out = y[global_var_ids[2]];
  double q_in = y[global_var_ids[1]];

  // Get parameters
  double Rmin = parameters[global_param_ids[ParamId::RMIN]];
  double Rmax = parameters[global_param_ids[ParamId::RMAX]];
  double epsilon = parameters[global_param_ids[ParamId::EPSILON]];

  // Nonlinear terms
  system.C(global_eqn_ids[0]) =
      -calculate_beta(p_in, p_out, epsilon, Rmax, Rmin);
  system.C(global_eqn_ids[2]) =
      -set_valve_status(p_in, p_out, epsilon, Rmax, Rmin);

  // Derivatives of non-linear terms
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[0]) =
      -calculate_dbeta(p_in, p_out, epsilon, Rmax, Rmin);
  system.dC_dy.coeffRef(global_eqn_ids[0], global_var_ids[2]) =
      calculate_dbeta(p_in, p_out, epsilon, Rmax, Rmin);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[0]) =
      -dvalve_status(p_in, p_out, epsilon, Rmax, Rmin);
  system.dC_dy.coeffRef(global_eqn_ids[2], global_var_ids[2]) =
      dvalve_status(p_in, p_out, epsilon, Rmax, Rmin);
}

double ValveSpline::calculate_beta(double p_in, double p_out, double epsilon,
                                   double R_max, double R_min) {
  double s = (p_in - p_out + 0.5 * epsilon) / epsilon;
  double p_0 = -0.5 * epsilon / R_max;
  double p_1 = 0.5 * epsilon / R_min;
  double m_0 = 1.0 / R_max;
  double m_1 = 1.0 / R_min;
  double h_00 = 2.0 * pow(s, 3) - 3.0 * pow(s, 2) + 1.0;
  double h_10 = pow(s, 3) - 2.0 * pow(s, 2) + s;
  double h_01 = -2.0 * pow(s, 3) + 3.0 * pow(s, 2);
  double h_11 = pow(s, 3) - pow(s, 2);

  if (p_in < p_out - 0.5 * epsilon) {
    return (p_in - p_out) / R_max;
  } else if (p_in < p_out + 0.5 * epsilon) {
    return h_00 * p_0 + h_10 * m_0 * epsilon + h_01 * p_1 +
           h_11 * m_1 * epsilon;
  } else {
    return (p_in - p_out) / R_min;
  }
}

double ValveSpline::calculate_dbeta(double p_in, double p_out, double epsilon,
                                    double R_max, double R_min) {
  double s = (p_in - p_out + 0.5 * epsilon) / epsilon;
  double p_0 = -0.5 * epsilon / R_max;
  double p_1 = 0.5 * epsilon / R_min;
  double m_0 = 1.0 / R_max;
  double m_1 = 1.0 / R_min;
  double h_00 = 2.0 * pow(s, 3) - 3.0 * pow(s, 2) + 1.0;
  double h_10 = pow(s, 3) - 2.0 * pow(s, 2) + s;
  double h_01 = -2.0 * pow(s, 3) + 3.0 * pow(s, 2);
  double h_11 = pow(s, 3) - pow(s, 2);

  if (p_in < p_out - 0.5 * epsilon) {
    return 1.0 / R_max;
  } else if (p_in < p_out + 0.5 * epsilon) {
    return (1 - s) * (1 / R_max) + s * (1 / R_min);
  } else {
    return 1.0 / R_min;
  }
}

double ValveSpline::set_valve_status(double p_in, double p_out, double epsilon,
                                     double R_max, double R_min) {
  double s = (p_in - p_out + 0.5 * epsilon) / epsilon;
  if (p_in < p_out - 0.5 * epsilon) {
    return 0;
  } else if (p_in < p_out + 0.5 * epsilon) {
    return ((1 - s) * (1 / R_max) + s * (1 / R_min) - 1 / R_max) /
           (1 / R_min - 1 / R_max);
  } else {
    return 1.0;
  }
}

double ValveSpline::dvalve_status(double p_in, double p_out, double epsilon,
                                  double R_max, double R_min) {
  double s = (p_in - p_out + 0.5 * epsilon) / epsilon;
  if (p_in < p_out - 0.5 * epsilon) {
    return 0;
  } else if (p_in < p_out + 0.5 * epsilon) {
    return 1 / epsilon;
  } else {
    return 0;
  }
}
