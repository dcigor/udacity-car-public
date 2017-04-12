#include "ukf.h"
#include "tools.h"
#include <iostream>
#include <cmath>
#include <initializer_list>

namespace {
    double clip_angle(double d) {
        while (d> M_PI) d-=2*M_PI;
        while (d<-M_PI) d+=2*M_PI;
        return d;
    }

    VectorXd clip_angle(const VectorXd &x, const std::initializer_list<size_t> &indices) {
        VectorXd d = x;
        for (size_t index : indices) {
            d(index) = clip_angle(d(index));
        }
        return d;
    }
}

// Initializes Unscented Kalman filter
UKF::UKF() {
    //set weights
    weights_[0] = lambda_aug_/(lambda_aug_+n_aug_);
    for (int i=1; i<2*n_aug_+1; ++i) {
        weights_[i] = .5/(lambda_aug_+n_aug_);
    }

    R_radar_(0,0) = std_radr_  *std_radr_;
    R_radar_(1,1) = std_radphi_*std_radphi_;
    R_radar_(2,2) = std_radrd_ *std_radrd_;

    //measurement covariance matrices
    R_laser_ << std_laspx_*std_laspx_, 0,
                0, std_laspy_*std_laspy_;

    //laser measurement matrix
    H_laser_ << 1, 0, 0, 0, 0,
                0, 1, 0, 0, 0;
    Ht_laser_ = H_laser_.transpose();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage &meas_package) {
    if (meas_package.sensor_type_ != MeasurementPackage::RADAR) {
        //return;
    }

    if (!is_initialized_) {
        init(meas_package);
        return;
    }

    // for small time increments: skip prediction step
    const double dt = std::max(0.00001, (meas_package.timestamp_ - previous_timestamp_) / 1000000.0);
    //if (dt > 0.00001)
    {
        Prediction(dt);
        previous_timestamp_ = meas_package.timestamp_;
    }

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    }   else {
        UpdateLidar(meas_package);
    }
}

void UKF::init(const MeasurementPackage &meas_package) {
    is_initialized_     = true;
    previous_timestamp_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Convert radar from polar to cartesian coordinates and initialize state.
        const double rho     = meas_package.raw_measurements_[0];
        const double phi     = clip_angle(meas_package.raw_measurements_[1]);
        const double rho_dot = meas_package.raw_measurements_[2];
        const double px      = cos(phi)*rho, py = sin(phi)*rho;
        is_phi_initialized_  = rho > 0.0001; // for zero position we have no initial info about the angle
        x_ << px, py, rho_dot, 0, 0; //[pos1 pos2 vel_abs yaw_angle yaw_rate]
    }   else {
        // Initialize state: initial position and zero speed
        //[pos1 pos2 vel_abs yaw_angle yaw_rate]
        x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(const double dt) {
    Xsig_pred_ = PredictSigmaPoints(dt);
    x_         = clip_angle(WeightedMean(Xsig_pred_), {3,4});
    P_         = PredictP(x_, Xsig_pred_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage &meas_package) {
    //update the state by using Kalman Filter equations
    const VectorXd z  = meas_package.raw_measurements_;
    const VectorXd y  = z - H_laser_ * x_;
    const MatrixXd S  = H_laser_ * P_ * Ht_laser_ + R_laser_;
    const MatrixXd Si = S.inverse();
    const MatrixXd K  = P_ * Ht_laser_ * Si;

    x_        += K*y;
    P_         = (I_-K*H_laser_)*P_;
    NIS_laser_ = y.transpose() * Si * y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage &meas_package) {
    const MatrixXd Zsig   = PredictRadarMeasurementSigmaPoints(Xsig_pred_);
    const VectorXd z_pred = clip_angle(WeightedMean(Zsig), {1});
          VectorXd y      = clip_angle(meas_package.raw_measurements_-z_pred, {1});
    const MatrixXd S      = RadarMeasurementCovarianceS(z_pred, Zsig);
    const MatrixXd Si     = S.inverse();
    const MatrixXd Tc     = RadarTc(x_, Xsig_pred_, z_pred, Zsig);
    const MatrixXd K      = Tc*Si;

    if (!is_phi_initialized_) {
        is_phi_initialized_ = true;
        x_(1) = y(1);
        y(1)  = 0;
    }

    //update state mean and covariance matrix
    std::cout << "UR S\n" << S << "\nSi\n" << Si << "\nTc\n" << Tc << "\nK\n" << K << "\n\n";
    std::cout << "UR x " << x_ << "\ny " << y << "\nkg " << (K*y) << "\n\n";

    x_         = clip_angle(x_+K*y, {3, 4});
    P_        -= K*S*K.transpose();
    NIS_radar_ = y.transpose() * Si * y;
}

MatrixXd UKF::AugmentedSigmaPoints() {
    //create augmented mean state
    VectorXd x_aug = VectorXd::Zero(7);
    x_aug.topLeftCorner(n_x_, 1) = x_;

    //create augmented covariance matrix
    MatrixXd P_aug = MatrixXd::Zero(7, 7);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_x_  , n_x_  ) = std_a_    *std_a_;
    P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;

    //create square root matrix
    const MatrixXd A = P_aug.llt().matrixL();

    //create augmented sigma points
    const double k    = sqrt(lambda_aug_+n_aug_);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
    Xsig_aug.col(0)   = x_aug;
    for (int i=0; i<n_aug_; ++i) {
        Xsig_aug.col(i+1       ) = clip_angle(x_aug+k*A.col(i), {3,4});
        Xsig_aug.col(i+1+n_aug_) = clip_angle(x_aug-k*A.col(i), {3,4});
    }

    return Xsig_aug;
}

// one column - one point.
MatrixXd UKF::PredictSigmaPoints(const double dt) {
    const MatrixXd Xsig_aug = AugmentedSigmaPoints();

    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd::Zero(n_x_, 2*n_aug_+1);
    for (int c = 0; c < 2*n_aug_+1; ++c) {
        const VectorXd &x = Xsig_aug.col(c); // one augmented point. 7 rows
        VectorXd     y    = x.head(n_x_);    // one predicted point. 5 rows
        const double v    = x(2);
        const double yaw  = clip_angle(x(3));
        const double yawd = clip_angle(x(4));
        const double nu   = x(5);
        const double nu2  = x(6);
        if (fabs(yawd) < 0.00001) {
            // straight motion
            y(0) += v*dt*cos(yaw);
            y(1) += v*dt*sin(yaw);
        }   else {
            // turning
            y(0) += v/yawd*( sin(yaw+yawd*dt)-sin(yaw)) + .5*dt*dt*cos(yaw)*nu;
            y(1) += v/yawd*(-cos(yaw+yawd*dt)+cos(yaw)) + .5*dt*dt*sin(yaw)*nu;
        }
        y(2) += dt*nu;
        y(3) += yawd*dt + 0.5*dt*dt*nu2;
        y(4) += dt*nu2;

        Xsig_pred.col(c) = clip_angle(y, {3,4});
    }

    std::cout << "PSP Xsig_aug\n" << Xsig_aug << std::endl << "Xsig_pred\n" << Xsig_pred << std::endl << std::endl;
    return Xsig_pred;
}

VectorXd UKF::WeightedMean(const MatrixXd &sig_points) {
    VectorXd x = weights_(0) * sig_points.col(0);
    for (int i = 1; i < 2*n_aug_+1; ++i) {
        x += weights_(i) * sig_points.col(i);
    }
    return x;
}

MatrixXd UKF::PredictP(const VectorXd &x, const MatrixXd &Xsig_pred) {
    //create covariance matrix for prediction
    MatrixXd P = MatrixXd::Zero(n_x_, n_x_);
    for (int i = 0; i < 2*n_aug_+1; ++i) {  //iterate over sigma points
        // state difference
        const VectorXd y = clip_angle(Xsig_pred.col(i)-x, {3});
        P += weights_(i) * y * y.transpose();
    }
    return P;
}

MatrixXd UKF::PredictRadarMeasurementSigmaPoints(const MatrixXd &Xsig_pred) {
    //transform sigma points into 3-d radar measurement space (r, phi, and r_dot)
    MatrixXd Zsig = MatrixXd::Zero(3, 2*n_aug_+1);

    for (int i=0; i < 2*n_aug_+1; ++i) {
        const double px  = Xsig_pred(0,i);
        const double py  = Xsig_pred(1,i);
        const double v   = Xsig_pred(2,i);
        const double k   = Xsig_pred(3,i);
        const double rho = sqrt(px*px+py*py);
        Zsig(0,i)        = rho;
        Zsig(1,i)        = atan2(py, px);
        Zsig(2,i)        = rho>0.000001 ? (px*cos(k)*v+py*sin(k)*v)/rho : 0;
    }
    return Zsig;
}

MatrixXd UKF::RadarMeasurementCovarianceS(const VectorXd &z_pred, const MatrixXd &Zsig) {
    MatrixXd S = R_radar_;
    for (int i=0; i < 2*n_aug_+1; ++i) {
        const VectorXd v = clip_angle(Zsig.col(i)-z_pred, {1});
        S += weights_(i)*v*v.transpose();
    }
    return S;
}

//calculate cross correlation matrix
MatrixXd UKF::RadarTc(const VectorXd &x, const MatrixXd &Xsig_pred, const VectorXd &z_pred, const MatrixXd &Zsig) {
    MatrixXd Tc = MatrixXd::Zero(n_x_, 3); //radar can measure 3 dimensions: r, phi, and r_dot
    for (int i=0; i<2*n_aug_+1; ++i) {
        Tc += weights_(i) * clip_angle(Xsig_pred.col(i)-x, {1}) * clip_angle(Zsig.col(i)-z_pred, {1}).transpose();
    }
    return Tc;
}
