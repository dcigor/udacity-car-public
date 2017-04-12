#include "ukf.h"
#include "tools.h"

// Initializes Unscented Kalman filter
UKF::UKF() {
    //set weights
    weights_[0] = lambda_aug_/(lambda_aug_+n_aug_);
    for (int i=1; i<2*n_aug_+1; ++i) {
        weights_[i] = .5/(lambda_aug_+n_aug_);
    }

    radar_R_(0,0) = std_radr_  *std_radr_;
    radar_R_(1,1) = std_radphi_*std_radphi_;
    radar_R_(2,2) = std_radrd_ *std_radrd_;
    
    //laser measurement matrix
    H_laser_ << 1, 0, 0, 0, 0,
                0, 1, 0, 0, 0;
    Ht_laser_ = H_laser_.transpose();
    
    //measurement covariance matrix - laser
    R_laser_ << std_laspx_, 0,
                0, std_laspy_;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage &meas_package) {
    if (!is_initialized_) {
        init(meas_package);
        return;
    }

    // for small time increments: skip prediction step
    const double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
    if (dt > 0.00001) {
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
        const double phi     = meas_package.raw_measurements_[1];
        const double rho_dot = meas_package.raw_measurements_[2];
        const double px      = cos(phi)*rho, py = sin(phi)*rho;
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
    x_         = WeightedMean(Xsig_pred_);
    P_         = PredictP(x_, Xsig_pred_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage &meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    //update the state by using Kalman Filter equations
    const VectorXd z = meas_package.raw_measurements_;
    const VectorXd y = z - H_laser_ * x_;
    const MatrixXd S = H_laser_ * P_ * Ht_laser_ + R_laser_;
    const MatrixXd K = P_ * Ht_laser_ * S.inverse();

    x_ += K*y;
    P_ = (I_-K*H_laser_)*P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage &meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    const MatrixXd Zsig   = PredictRadarMeasurementSigmaPoints(Xsig_pred_);
    const VectorXd z_pred = WeightedMean(Zsig);
    const MatrixXd S      = RadarMeasurementCovarianceS(z_pred, Zsig);
    const MatrixXd Tc     = RadarTc(x_, Xsig_pred_, z_pred, Zsig);
    
    //calculate Kalman gain K;
    const MatrixXd K = Tc*S.inverse();
    
    //update state mean and covariance matrix
    const VectorXd y = meas_package.raw_measurements_-z_pred;
    x_ += K*y;
    P_ -= K*S*K.transpose();
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
        Xsig_aug.col(i+1       ) = x_aug+k*A.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug-k*A.col(i);
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
        const double yaw  = x(3);
        const double yawd = x(4);
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

        Xsig_pred.col(c) = y;
    }

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
        VectorXd x_diff = Xsig_pred.col(i) - x;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P += weights_(i) * x_diff * x_diff.transpose();
    }
    return P;
}

MatrixXd UKF::PredictRadarMeasurementSigmaPoints(const MatrixXd &Xsig_pred) {
    //transform sigma points into 3-d radar measurement space (r, phi, and r_dot)
    MatrixXd Zsig = MatrixXd::Zero(3, 2*n_aug_+1);

    for (int i=0; i < 2*n_aug_+1; ++i) {
        const double px = Xsig_pred(0,i);
        const double py = Xsig_pred(1,i);
        const double v  = Xsig_pred(2,i);
        const double k  = Xsig_pred(3,i);
        Zsig(0,i) = sqrt(px*px+py*py);
        Zsig(1,i) = atan2(py, px);
        Zsig(2,i) = (px*cos(k)*v+py*sin(k)*v)/Zsig(0,i);
    }
    return Zsig;
}

MatrixXd UKF::RadarMeasurementCovarianceS(const VectorXd &z_pred, const MatrixXd &Zsig) {
    MatrixXd S = radar_R_;
    for (int i=0; i < 2*n_aug_+1; ++i) {
        const VectorXd v = Zsig.col(i)-z_pred;
        S += weights_(i)*v*v.transpose();
    }
    return S;
}

//calculate cross correlation matrix
MatrixXd UKF::RadarTc(const VectorXd &x, const MatrixXd &Xsig_pred, const VectorXd &z_pred, const MatrixXd &Zsig) {
    MatrixXd Tc = MatrixXd::Zero(n_x_, 3); //radar can measure 3 dimensions: r, phi, and r_dot
    for (int i=0; i<2*n_aug_+1; ++i) {
        Tc += weights_(i)*(Xsig_pred.col(i)-x)*(Zsig.col(i)-z_pred).transpose();
    }
    return Tc;
}
