#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Geometry"

using CppAD::AD;

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

namespace {
    
    // number of time steps. for values of 25 and higher the solver cannot keep up.
    // smaller values do not let the car to see the turns, until its too late.
    const size_t N = 18;

    // for values 0.07 and higher - the solver cannot keep up. apparently it affects solver's time.
    // smaller values do not let the car to see the turns.
    const double dt    = 0.05; // time step duration
    
    // the highest value allowed by the simulator.
    // it might be possible to achieve higher speed after testing in updated simulator.
    const double ref_v = 100;  // The reference velocity, in mph.
    
    // This value assumes the model presented in the classroom is used.
    //
    // It was obtained by measuring the radius formed by running the vehicle in the
    // simulator around in a circle with a constant steering angle and velocity on a
    // flat terrain.
    //
    // Lf was tuned until the the radius formed by the simulating the model
    // presented in the classroom matched the previous radius.
    //
    // This is the length from front to CoG that has a similar radius.
    const double Lf = 2.67;

    // The solver takes all the state variables and actuator
    // variables in a singular vector. Thus, we should to establish
    // when one variable starts and another ends to make our lifes easier.
    const size_t x_start     = 0;
    const size_t y_start     = x_start     + N;
    const size_t psi_start   = y_start     + N;
    const size_t v_start     = psi_start   + N;
    const size_t cte_start   = v_start     + N;
    const size_t epsi_start  = cte_start   + N;
    const size_t delta_start = epsi_start  + N;
    const size_t a_start     = delta_start + N - 1;
}

// standard rigid body 2d transformation: translation and rotation
World::World(const Eigen::VectorXd &state) {
    const double x   = state(0);
    const double y   = state(1);
    const double yaw = state(2);

    const auto rot = Eigen::Rotation2Dd(yaw);
    const auto R   = rot.matrix();
    const auto Rt  = R.transpose();
    const Eigen::Vector2d translation(x,y);
    
    X_.topLeftCorner (2,2) = R;
    X_.topRightCorner(2,1) = translation;
    
    X_inv_.topLeftCorner (2,2) = Rt;
    X_inv_.topRightCorner(2,1) = -Rt*translation;
}

// `fg` is a vector containing the cost and constraints.
// `vars` is a vector containing the variable values (state & actuators).
void FG_eval::operator()(ADvector& fg, const ADvector& vars) {
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // The part of the cost based on the reference state.
    for (int t = 0; t < N; t++) {
        fg[0] += 500*CppAD::pow(vars[cte_start  + t], 2);
        fg[0] += 200*CppAD::pow(vars[epsi_start + t], 2);
        fg[0] +=     CppAD::pow(vars[v_start    + t] - ref_v, 2);

        // reduce centrifugal force: v*steering:
        // avoid turning at high speeds.
        fg[0] += 30*CppAD::pow(vars[epsi_start+t]*vars[v_start+t], 2);
    }
    
    // Minimize the use of actuators.
    for (int t = 0; t < N - 1; t++) {
        fg[0] += 100*CppAD::pow(vars[delta_start + t], 2);
        fg[0] +=     CppAD::pow(vars[a_start     + t], 2);
    }
    
    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < N - 2; t++) {
        fg[0] += 400*CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
        fg[0] +=     CppAD::pow(vars[a_start     + t + 1] - vars[a_start     + t], 2);
    }
    
    // Setup Constraints
    
    // Initial constraints
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start   ] = vars[x_start   ];
    fg[1 + y_start   ] = vars[y_start   ];
    fg[1 + psi_start ] = vars[psi_start ];
    fg[1 + v_start   ] = vars[v_start   ];
    fg[1 + cte_start ] = vars[cte_start ];
    fg[1 + epsi_start] = vars[epsi_start];
    
    // The rest of the constraints
    for (size_t t = 1; t < N; t++) {
        // The state at time t+1 .
        AD<double> x1    = vars[x_start    + t];
        AD<double> y1    = vars[y_start    + t];
        AD<double> psi1  = vars[psi_start  + t];
        AD<double> v1    = vars[v_start    + t];
        AD<double> cte1  = vars[cte_start  + t];
        AD<double> epsi1 = vars[epsi_start + t];
        
        // The state at time t.
        AD<double> x0    = vars[x_start    + t - 1];
        AD<double> y0    = vars[y_start    + t - 1];
        AD<double> psi0  = vars[psi_start  + t - 1];
        AD<double> v0    = vars[v_start    + t - 1];
        AD<double> cte0  = vars[cte_start  + t - 1];
        AD<double> epsi0 = vars[epsi_start + t - 1];
        
        // Only consider the actuation at time t.
        AD<double> delta0 = vars[delta_start + t - 1];
        AD<double> a0     = vars[a_start     + t - 1];
        
        /// we approximate using the cubic polyline. here is the current value: y = a*x^3 + b*x^2 + c*x + d
        AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
        /// here is the first derivative (tangent): y' = 3*a*x^2 + 2*b*x + c
        AD<double> psides0 = CppAD::atan(coeffs[1] + 2*coeffs[2]*x0 + 3*coeffs[3]*x0*x0);

        // Here's `x` to get you started.
        // The idea here is to constraint this value to be 0.
        //
        // Recall the equations for the model:
        // x_  [t+1] = x[t] + v[t] * cos(psi[t]) * dt
        // y_  [t+1] = y[t] + v[t] * sin(psi[t]) * dt
        // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
        // v_  [t+1] = v[t] + a[t] * dt
        // cte [t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
        // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
        fg[1 + x_start    + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
        fg[1 + y_start    + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
        fg[1 + psi_start  + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
        fg[1 + v_start    + t] = v1 - (v0 + a0 * dt);
        fg[1 + cte_start  + t] = (cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt)));
        fg[1 + epsi_start + t] = (epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt));
    }
}

// MPC class definition
void MPC::Solve(const double x, const double y, const double psi, const double v, const double cte, const double epsi, const Eigen::VectorXd coeffs) {
    // number of independent variables: N timesteps == N - 1 actuations
    const size_t n_vars = N * 6 + (N - 1) * 2;
    const size_t n_constraints = N * 6; // Number of constraints

    // Initial value of the independent variables.
    // Should be 0 except for the initial values.
    Dvector vars(n_vars);
    for (size_t i = 0; i < n_vars; i++) {
        vars[i] = 0.0;
    }
    // Set the initial variable values
    vars[x_start]    = x;
    vars[y_start]    = y;
    vars[psi_start]  = psi;
    vars[v_start]    = v;
    vars[cte_start]  = cte;
    vars[epsi_start] = epsi;

    // Lower and upper limits for x
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);
    
    // Set all non-actuators upper and lower limits
    // to the max negative and positive values.
    for (size_t i = 0; i < delta_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    // The upper and lower limits of delta
    for (size_t i = delta_start; i < a_start; i++) {
        vars_lowerbound[i] = -M_PI/10;
        vars_upperbound[i] =  M_PI/10;
    }
    
    // Acceleration/decceleration upper and lower limits.
    for (size_t i = a_start; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }
    
    // Lower and upper limits for constraints
    // All of these should be 0 except the initial state indices.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (size_t i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
    constraints_lowerbound[x_start   ] = x;
    constraints_lowerbound[y_start   ] = y;
    constraints_lowerbound[psi_start ] = psi;
    constraints_lowerbound[v_start   ] = v;
    constraints_lowerbound[cte_start ] = cte;
    constraints_lowerbound[epsi_start] = epsi;

    constraints_upperbound[x_start   ] = x;
    constraints_upperbound[y_start   ] = y;
    constraints_upperbound[psi_start ] = psi;
    constraints_upperbound[v_start   ] = v;
    constraints_upperbound[cte_start ] = cte;
    constraints_upperbound[epsi_start] = epsi;

    // Object that computes objective and constraints
    FG_eval fg_eval(coeffs);

    // options
    std::string options;
    options += "Integer print_level  0\n";
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(options, vars, vars_lowerbound, vars_upperbound,
        constraints_lowerbound, constraints_upperbound, fg_eval, solution_);

    // Check some of the solution values
    const bool ok = solution_.status == CppAD::ipopt::solve_result<Dvector>::success;
}

vector<double> MPC::getRow(const size_t row) const {
    vector<double> res;
    const size_t n_vars = N * 6 + (N - 1) * 2;
    for (size_t i=0; i<N; ++i) {
        res.push_back(solution_.x[row*N+i]);
    }
    return res;
}

double MPC::getSteeringDelta() const {
    return solution_.x[delta_start];
}

double MPC::getAcceleration() const {
    return solution_.x[a_start];
}

double MPC::getCost() const {
    return solution_.obj_value;
}
