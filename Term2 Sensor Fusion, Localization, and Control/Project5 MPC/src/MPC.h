#ifndef MPC_H
#define MPC_H

#include <vector>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

// world transformation
// local coordinates: coordinates from iside the car. x - points forward, y - points left to the car
class World {
public:
    World(double x, double y, double yaw);
    
    // convert global coordinates to local (car) coordinates
    Eigen::Vector3d globalToLocal(const Eigen::Vector3d &v) const {
        return X_inv_ * v;
    }
    Eigen::Vector3d globalToLocal(const double x, const double y) const {
        return globalToLocal(Eigen::Vector3d(x,y,1));
    }
    vector<Eigen::Vector3d> globalToLocal(const vector<double> &x, const vector<double> &y) const {
        vector<Eigen::Vector3d> res(x.size());
        for (size_t i=0; i<x.size(); ++i) {
            res[i] = globalToLocal(x[i], y[i]);
        }
        return res;
    }

    // convert local coordinates to global
    Eigen::Vector3d localToGlobal(const Eigen::Vector3d &v) const {
        return X_ * v;
    }
    Eigen::Vector3d localToGlobal(const double x, const double y) const {
        return localToGlobal(Eigen::Vector3d(x,y,1));
    }

    // utility method that should be moved sowhere else.
    // get a vector of x-coordinates or y-coordinates from a vector of vectors
    static vector<double> getRow(const size_t row, const vector<Eigen::Vector3d> &v) {
        vector<double> res(v.size());
        for (size_t i=0; i<v.size(); ++i) {
            res[i] = v[i](row);
        }
        return res;
    }
private:
    // world transform matrix and its inverse
    Eigen::Matrix3d X_     = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d X_inv_ = Eigen::Matrix3d::Identity();
};

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
    void Solve(double x, double y, double psi, double v, double cte, double epsi, Eigen::VectorXd coeffs);

    // accessors
    vector<double> getRow(size_t row) const;
    double getSteeringDelta() const;
    double getAcceleration () const;
    double getCost         () const;
private:
    typedef CPPAD_TESTVECTOR(double) Dvector;
    CppAD::ipopt::solve_result<Dvector> solution_;
};

#endif /* MPC_H */
