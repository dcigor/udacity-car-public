#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
public:
    KalmanFilter(float process_noise);

    void initX(float px, float py, float vx, float vy);
    Eigen::VectorXd x() const {return x_;}
    Eigen::MatrixXd P() const {return P_;}

    /**
     * Prediction Predicts the state and the state covariance
     * using the process model
     * @param dt Time between k and k+1 in s
     */
    void Predict(float dt, float noise_ax, float noise_ay);

    /**
     * Updates the state by using standard Kalman Filter equations
     * @param R Measurement covariance matrix
     * @param z The measurement at k+1
     */
    void Update(const Eigen::MatrixXd &R, const Eigen::VectorXd &z);

    /**
     * Updates the state by using Extended Kalman Filter equations
     * @param R Measurement covariance matrix
     * @param z The measurement at k+1
     */
    void UpdateEKF(const Eigen::MatrixXd &R, const Eigen::VectorXd &z);
private:
    // state vector
    Eigen::VectorXd x_ = Eigen::VectorXd(4);

    // state covariance matrix
    Eigen::MatrixXd P_ = Eigen::MatrixXd(4, 4);

    // state transistion matrix
    Eigen::MatrixXd F_ = Eigen::MatrixXd::Identity(4, 4);

    // process covariance matrix
    Eigen::MatrixXd Q_ = Eigen::MatrixXd(4, 4);

    // measurement matrix
    Eigen::MatrixXd H_ = Eigen::MatrixXd(2, 4), Ht_;

    const Eigen::MatrixXd I_ = Eigen::MatrixXd::Identity(4, 4);
};

#endif /* KALMAN_FILTER_H_ */
