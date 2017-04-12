#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
    UKF();

    // public accessors to the members
    const VectorXd &x() const {return x_;}
    double NIS_laser () const {return NIS_laser_;}
    double NIS_radar () const {return NIS_radar_;}
    
    /**
     * ProcessMeasurement
     * @param meas_package The latest measurement data of either radar or laser
     */
    void ProcessMeasurement(const MeasurementPackage &meas_package);
private:
    void init(const MeasurementPackage &meas_package);
    
    /**
     * Prediction Predicts sigma points, the state, and the state covariance
     * matrix
     * @param delta_t Time between k and k+1 in s
     */
    void Prediction(double delta_t);

    /**
     * Updates the state and the state covariance matrix using a laser measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateLidar(const MeasurementPackage &meas_package);

    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateRadar(const MeasurementPackage &meas_package);
    
    MatrixXd AugmentedSigmaPoints();
    MatrixXd PredictSigmaPoints(double dt);
    VectorXd WeightedMean(const MatrixXd &sig_points);
    MatrixXd PredictP(const VectorXd &x, const MatrixXd &Xsig_pred);

    MatrixXd PredictRadarMeasurementSigmaPoints(const MatrixXd &Xsig_pred);
    MatrixXd RadarMeasurementCovarianceS(const VectorXd &z_pred, const MatrixXd &Zsig);
    MatrixXd RadarTc(const VectorXd &x, const MatrixXd &Xsig_pred, const VectorXd &z_pred, const MatrixXd &Zsig);

    // initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_ = false, is_phi_initialized_ = false;
    
    // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_ = VectorXd::Zero(5);
    
    // state covariance matrix
    MatrixXd P_ = MatrixXd::Identity(5, 5);

    // predicted sigma points matrix
    MatrixXd Xsig_pred_;

    // State dimension & Augmented state dimension
    const int n_x_ = 5, n_aug_ = 7;

    MatrixXd R_radar_ = MatrixXd::Zero(3,3);

    // Weights of sigma points
    VectorXd weights_ = VectorXd::Zero(2*n_aug_+1);

    // Sigma point spreading parameter
    const double lambda_ = 3-n_x_, lambda_aug_ = 3-n_aug_;

    // the current NIS for radar and laser
    double NIS_radar_ = 0, NIS_laser_ = 0;

    long previous_timestamp_ = 0;

    // process noise: needs tuning
    // Process noise standard deviation longitudinal acceleration in m/s^2
    const double std_a_ = 0.5;
    // Process noise standard deviation yaw acceleration in rad/s^2
    const double std_yawdd_ = 2;

    // measurement noise: provided by the device manufacturer. do not change
    // Laser measurement noise standard deviation px and py in m
    const double std_laspx_ = 0.15, std_laspy_ = std_laspx_;
    // Radar measurement noise standard deviation radius in m
    const double std_radr_ = 0.3;
    // Radar measurement noise standard deviation angle in rad
    const double std_radphi_ = 0.03;
    // Radar measurement noise standard deviation radius change in m/s
    const double std_radrd_ = 0.3;

    // laser
    const MatrixXd I_ = MatrixXd::Identity(5, 5);
    MatrixXd H_laser_ = MatrixXd(2, 5), Ht_laser_; // measurement matrix
    MatrixXd R_laser_ = MatrixXd(2, 2); // measurement covariance
};

#endif /* UKF_H */
