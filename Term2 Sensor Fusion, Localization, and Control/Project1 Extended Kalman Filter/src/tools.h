#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

namespace Tools {
    // A helper method to calculate RMSE.
     Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

    // A helper method to calculate Jacobians.
    std::pair<Eigen::MatrixXd, bool> CalculateJacobian(const Eigen::VectorXd& x);
}

#endif /* TOOLS_H_ */
