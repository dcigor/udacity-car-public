#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
:   ekf_(1000)
{
    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    const bool isRadar = measurement_pack.sensor_type_ == MeasurementPackage::RADAR;
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        // first measurement
        if (isRadar) {
            // Convert radar from polar to cartesian coordinates and initialize state.
            const float rho     = measurement_pack.raw_measurements_[0];
            const float phi     = measurement_pack.raw_measurements_[1];
            const float rho_dot = measurement_pack.raw_measurements_[2];
            const float px      = cos(phi)*rho,     py = sin(phi)*rho;
            const float vx      = cos(phi)*rho_dot, vy = sin(phi)*rho_dot;
            ekf_.initX(px, py, vx, vy);
        }   else {
            // Initialize state: initial position and zero speed
            ekf_.initX(measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0);
        }

        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     */
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    if (dt < 0.000001) {
        dt = 0.000001; // fallback to the smallest possible time increment
    }

    ekf_.Predict(dt, 9, 9); // Use noise_ax = 9 and noise_ay = 9 for your Q matrix.

    /*****************************************************************************
     *  Update
     ****************************************************************************/
    if (isRadar) {
        ekf_.UpdateEKF(R_radar_, measurement_pack.raw_measurements_);
    }   else {
        ekf_.Update(R_laser_, measurement_pack.raw_measurements_);
    }

    previous_timestamp_ = measurement_pack.timestamp_;
}
