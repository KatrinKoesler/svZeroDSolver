// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
/**
 * @file ValveSpline.h
 * @brief model::ValveSpline source file
 */
#ifndef SVZERODSOLVER_MODEL_VALVESPLINE_HPP_
#define SVZERODSOLVER_MODEL_VALVESPLINE_HPP_

#include <math.h>

#include "Block.h"
#include "SparseSystem.h"
#include "debug.h"

/**
 * @brief Valve (spline) block.
 *
 * Models the pressure-flow relationship across a valve using a cubic Hermite
 * spline to smoothly interpolate between a closed state (resistance \f$R_{max}
 * \f$) and an open state (resistance \f$R_{min}\f$). The transition is governed
 * by the pressure difference \f$\Delta p = P_{in} - P_{out}\f$ over an interval
 * of width \f$\epsilon\f$ (Epsilon). See Hirschvogel et al. 2024.
 *
 * \f[
 * \begin{circuitikz} \draw
 * node[left] {$Q_{in}$} [-latex] (0,0) -- (0.8,0);
 * \draw (1,0) node[anchor=south]{$P_{in}$}
 * to [D, l=$\beta$, *-*] (3,0)
 * node[anchor=south]{$P_{out}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * Q_{in} - \beta(P_{in}, P_{out}) = 0
 * \f]
 *
 * \f[
 * Q_{in} - Q_{out} = 0
 * \f]
 *
 * \f[
 * \text{valve\_status} - \sigma(P_{in}, P_{out}) = 0
 * \f]
 *
 * where \f$\beta\f$ is the spline-interpolated flow and \f$\sigma \in [0,1]\f$
 * is the valve openness (0 = closed, 1 = open).
 *
 * Let \f$s = \frac{\Delta p + \epsilon/2}{\epsilon}\f$ with
 * \f$\Delta p = P_{in} - P_{out}\f$.
 *
 * \f[
 * \beta(\Delta p) =
 * \begin{cases}
 * \Delta p \,/\, R_{max} & \Delta p < -\epsilon/2 \\[4pt]
 * h_{00}(s)\,p_0 + h_{10}(s)\,m_0\,\epsilon
 *   + h_{01}(s)\,p_1 + h_{11}(s)\,m_1\,\epsilon
 *   & |\Delta p| \leq \epsilon/2 \\[4pt]
 * \Delta p \,/\, R_{min} & \Delta p > \epsilon/2
 * \end{cases}
 * \f]
 *
 * with cubic Hermite basis functions
 * \f$h_{00}=2s^3-3s^2+1\f$, \f$h_{10}=s^3-2s^2+s\f$,
 * \f$h_{01}=-2s^3+3s^2\f$, \f$h_{11}=s^3-s^2\f$,
 * boundary values \f$p_0 = -\epsilon/(2R_{max})\f$,
 * \f$p_1 = \epsilon/(2R_{min})\f$,
 * and slopes \f$m_0 = 1/R_{max}\f$, \f$m_1 = 1/R_{min}\f$.
 *
 * \f[
 * \sigma(\Delta p) =
 * \begin{cases}
 * 0 & \Delta p < -\epsilon/2 \\
 * s & |\Delta p| \leq \epsilon/2 \\
 * 1 & \Delta p > \epsilon/2
 * \end{cases}
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{lllll}P_{in} & Q_{in} &
 * P_{out} & Q_{out} & \text{valve\_status} \end{array}\right]^{T}
 * \f]
 *
 * \f[
 * \mathbf{E}^{e}=\mathbf{0}
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ccccc}
 * 0 & 1 & 0 &  0 & 0 \\
 * 0 & 1 & 0 & -1 & 0 \\
 * 0 & 0 & 0 &  0 & 1
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{c}^{e}=\left[\begin{array}{c}
 * -\beta(P_{in}, P_{out}) \\
 * 0 \\
 * -\sigma(P_{in}, P_{out})
 * \end{array}\right]
 * \f]
 *
 * \f[
 * \left(\frac{\partial\mathbf{c}}{\partial\mathbf{y}}\right)^{e} =
 * \left[\begin{array}{ccccc}
 * -\beta' & 0 &  \beta' & 0 & 0 \\
 *  0      & 0 &  0      & 0 & 0 \\
 * -\sigma' & 0 & \sigma' & 0 & 0
 * \end{array}\right]
 * \f]
 *
 * where the derivatives with respect to \f$\Delta p\f$ are:
 *
 * \f[
 * \beta'(\Delta p) =
 * \begin{cases}
 * 1/R_{max} & \Delta p < -\epsilon/2 \\
 * (1-s)/R_{max} + s/R_{min} & |\Delta p| \leq \epsilon/2 \\
 * 1/R_{min} & \Delta p > \epsilon/2
 * \end{cases}
 * \f]
 *
 * \f[
 * \sigma'(\Delta p) =
 * \begin{cases}
 * 0 & \Delta p < -\epsilon/2 \\
 * 1/\epsilon & |\Delta p| \leq \epsilon/2 \\
 * 0 & \Delta p > \epsilon/2
 * \end{cases}
 * \f]
 *
 * \f[
 * \left(\frac{\partial\mathbf{c}}{\partial\dot{\mathbf{y}}}\right)^{e} =
 * \mathbf{E}^{e}=\mathbf{0}
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block (parameters here from
 * Hirschvogel et al. 2024)
 *
 * * `0` Rmax: Maximum (closed) valve resistance
 * * `1` Rmin: Minimum (open) valve resistance
 * * `2` Epsilon: Pressure interval width for spline interpolation between
 *   closed and open valve states; calibrate to pressure differential expected
 * to be available across valve
 * * `3` upstream_block: Name of block connected upstream
 * * `4` downstream_block: Name of block connected downstream
 *
 * ### Usage in json configuration file
 *
 *     "valves": [
 *         {
 *             "type": "ValveSpline",
 *             "name": "valve",
 *             "params": {
 *                 "Rmax": 10.0e12,
 *                 "Rmin": 1.0e6,
 *                 "Epsilon": 0.5e3,
 *                 "upstream_block": "upstream_vessel",
 *                 "downstream_block": "downstream_vessel"
 *             }
 *         }
 *     ]
 *
 * ### Internal variables
 *
 * Name | Symbol | Unit | Description
 * -----|--------|------|------------
 * valve_status | \f$\sigma\f$ | - | Valve openness (0 = closed, 1 = fully open)
 *
 */
class ValveSpline : public Block {
 public:
  /**
   * @brief Local IDs of the parameters
   *
   */
  enum ParamId {
    RMAX = 0,
    RMIN = 1,
    EPSILON = 2,
  };

  /**
   * @brief Construct a new ValveSpline object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ValveSpline(int id, Model* model)
      : Block(id, model, BlockType::valve_spline, BlockClass::valve,
              {{"Rmax", InputParameter()},
               {"Rmin", InputParameter()},
               {"Epsilon", InputParameter()},
               {"upstream_block", InputParameter(false, false, false)},
               {"downstream_block", InputParameter(false, false, false)}}) {}

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set global_var_ids and global_eqn_ids of the element based on the
   * number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler& dofhandler) override;

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * Populates the solution-independent entries of \f$\mathbf{F}\f$.
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem& system,
                       std::vector<double>& parameters) override;

  /**
   * @brief Update the solution-dependent contributions of the element in a
   * sparse system
   *
   * Evaluates the nonlinear terms \f$\mathbf{c}\f$ and their Jacobian
   * \f$\partial\mathbf{c}/\partial\mathbf{y}\f$ at the current solution.
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Current solution
   * @param dy Current derivative of the solution
   */
  void update_solution(
      SparseSystem& system, std::vector<double>& parameters,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& dy) override;

  /**
   * @brief Compute the spline-interpolated flow \f$\beta(P_{in}, P_{out})\f$
   *
   * Uses a cubic Hermite spline to transition smoothly between the closed
   * (\f$R_{max}\f$) and open (\f$R_{min}\f$) resistance regimes over the
   * pressure interval of width \f$\epsilon\f$.
   *
   * @param p_in Inlet pressure
   * @param p_out Outlet pressure
   * @param epsilon Transition interval width
   * @param R_max Closed valve resistance
   * @param R_min Open valve resistance
   * @return double Flow \f$\beta\f$
   */
  double calculate_beta(double p_in, double p_out, double epsilon, double R_max,
                        double R_min);

  /**
   * @brief Compute \f$\beta'(\Delta p)\f$: derivative of \f$\beta\f$ with
   * respect to the pressure difference \f$\Delta p = P_{in} - P_{out}\f$
   *
   * In the transition zone this equals the linearly interpolated conductance
   * \f$(1-s)/R_{max} + s/R_{min}\f$.
   *
   * @param p_in Inlet pressure
   * @param p_out Outlet pressure
   * @param epsilon Transition interval width
   * @param R_max Closed valve resistance
   * @param R_min Open valve resistance
   * @return double \f$\mathrm{d}\beta/\mathrm{d}(\Delta p)\f$
   */
  double calculate_dbeta(double p_in, double p_out, double epsilon,
                         double R_max, double R_min);

  /**
   * @brief Compute the valve openness status \f$\sigma(P_{in}, P_{out})\f$
   *
   * Returns 0 when the valve is closed, 1 when fully open, and the
   * normalized spline parameter \f$s\f$ in the transition zone.
   *
   * @param p_in Inlet pressure
   * @param p_out Outlet pressure
   * @param epsilon Transition interval width
   * @param R_max Closed valve resistance
   * @param R_min Open valve resistance
   * @return double Valve status \f$\sigma \in [0, 1]\f$
   */
  double set_valve_status(double p_in, double p_out, double epsilon,
                          double R_max, double R_min);

  /**
   * @brief Compute \f$\sigma'(\Delta p)\f$: derivative of valve status with
   * respect to the pressure difference \f$\Delta p = P_{in} - P_{out}\f$
   *
   * Returns \f$1/\epsilon\f$ in the transition zone and 0 otherwise.
   *
   * @param p_in Inlet pressure
   * @param p_out Outlet pressure
   * @param epsilon Transition interval width
   * @param R_max Closed valve resistance
   * @param R_min Open valve resistance
   * @return double \f$\mathrm{d}\sigma/\mathrm{d}(\Delta p)\f$
   */
  double dvalve_status(double p_in, double p_out, double epsilon, double R_max,
                       double R_min);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{4, 0, 2};
};

#endif  // SVZERODSOLVER_MODEL_VALVESPLINE_HPP_
