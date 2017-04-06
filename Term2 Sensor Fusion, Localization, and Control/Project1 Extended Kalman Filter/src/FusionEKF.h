#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
    /**
     * Constructor.
     */
    FusionEKF();

    /**
     * Run the whole flow of the Kalman Filter from here.
     */
    void ProcessMeasurement(const MeasurementPackage &measurement_pack);

    /**
     * Kalman Filter update and prediction math lives in here.
     */
    KalmanFilter ekf_;

private:
    // check whether the tracking toolbox was initialized or not (first measurement)
    bool is_initialized_= false;

    // previous timestamp
    long previous_timestamp_ = 0;

    Eigen::MatrixXd R_laser_ = Eigen::MatrixXd(2, 2);
    Eigen::MatrixXd R_radar_ = Eigen::MatrixXd(3, 3);
    Eigen::MatrixXd H_laser_ = Eigen::MatrixXd(2, 4);
    Eigen::MatrixXd Hj_      = Eigen::MatrixXd(3, 4);
};

#endif /* FusionEKF_H_ */
