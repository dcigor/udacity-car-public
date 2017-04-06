#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse = VectorXd::Zero(4);

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    // ... your code here
    if (estimations.empty()) {
        std::cout << "CalculateRMSE() - Error - Division by Zero" << std::endl;
        return rmse;
    }

    //accumulate squared residuals
    for (size_t i=0; i<estimations.size(); ++i) {
        // ... your code here
        VectorXd d = estimations[i]-ground_truth[i];
        d = d.array()*d.array();
        rmse += d;
    }

    //calculate the mean
    rmse /= estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

std::pair<Eigen::MatrixXd, bool> Tools::CalculateJacobian(const VectorXd& x) {
    MatrixXd Hj(3,4);
    //recover state parameters
    const float px = x(0);
    const float py = x(1);
    const float vx = x(2);
    const float vy = x(3);

    //pre-compute a set of terms to avoid repeated calculation
    const float c1 = px*px+py*py;
    const float c2 = sqrt(c1);
    const float c3 = c1*c2;

    //check division by zero
    if (c1 < 0.00001) {
        std::cout << "CalculateJacobian() - Avoided Division by Zero" << std::endl;
        return {Hj, false};
    }

    //compute the Jacobian matrix
    Hj << px/c2, py/c2, 0, 0,
		  -py/c1, px/c1, 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    return {Hj, true};
}
