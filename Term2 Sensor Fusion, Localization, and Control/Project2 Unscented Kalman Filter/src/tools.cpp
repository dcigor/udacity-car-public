#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

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
