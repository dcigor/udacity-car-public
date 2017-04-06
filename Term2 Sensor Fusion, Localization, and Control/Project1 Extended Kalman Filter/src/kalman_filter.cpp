#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {
    template<typename T>
    T clip_angle(T angle) {
        while (angle < -M_PI) {
            angle += M_PI*2;
        }
        while (angle > M_PI) {
            angle -= M_PI*2;
        }
        return angle;
    }
    
    // function h converts position and speed from cartesian to polar
    VectorXd h(const VectorXd &x) {
        const float radius = sqrt(x[0]*x[0]+x[1]*x[1]);
        VectorXd res = VectorXd(3);
        res << radius,
               atan2(x[1], x[0]),
               (x[0]*x[2]+x[1]*x[3])/radius;
        return res;
    }
}

KalmanFilter::KalmanFilter(const float process_noise) {
    P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, process_noise, 0,
          0, 0, 0, process_noise;

    //measurement matrix
    H_ << 1, 0, 0, 0,
          0, 1, 0, 0;
    Ht_ = H_.transpose();
}

void KalmanFilter::initX(const float px, const float py, const float vx, const float vy) {
    x_ << px, py, vx, vy;
}

void KalmanFilter::Predict(const float dt, const float noise_ax, const float noise_ay) {
    const float dt_2 = dt   * dt;
    const float dt_3 = dt_2 * dt;
    const float dt_4 = dt_3 * dt;
    
    //Modify the F matrix so that the time is integrated
    F_(0,2) = dt;
    F_(1,3) = dt;
    
    // Set the process covariance matrix Q
    Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
          0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
          dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
          0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

    // Predict
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const Eigen::MatrixXd &R, const VectorXd &z) {
    //update the state by using Kalman Filter equations
    const VectorXd y = z - H_ * x_;
    const MatrixXd S = H_ * P_ * Ht_ + R;
    const MatrixXd K = P_ * Ht_ * S.inverse();

    x_ += K*y;
    P_ = (I_-K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const Eigen::MatrixXd &R, const VectorXd &z) {
    //update the state by using Extended Kalman Filter equations
    VectorXd y = z-h(x_);
    y[0] = clip_angle(y[0]);

    const auto Hj = Tools::CalculateJacobian(x_);
    if (!Hj.second) {
        return;
    }
    const MatrixXd Hjt = Hj.first.transpose();
    const MatrixXd S   = Hj.first * P_ * Hjt + R;
    const MatrixXd K   = P_ * Hjt * S.inverse();

    x_ += K*y;
    P_ = (I_-K*Hj.first)*P_;
}
