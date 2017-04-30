#include "PID.h"
#include <time.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::Init(double in_Kp, double in_Ki, double in_Kd) {
    Kp = in_Kp;
    Kd = in_Kd;
    Ki = in_Ki;
}

void PID::UpdateError(const double cte) {
    history_.push_back(cte);
    if (history_.size() > 100) {
        history_.erase(history_.begin());
    }

    p_error = cte;
    d_error = cte-prev_cte_;
    i_error = std::accumulate(history_.begin(), history_.end(), 0);

    prev_cte_ = cte;

    learn_coefficient(cte, Kp);
}

double PID::TotalError() {
    return -(p_error * Kp + d_error * Kd + i_error * Ki);
}

// this is only needed to find the best values of the coefficients.
void PID::learn_coefficient(const double cte, double &K) {
    const double steering = TotalError();
    const double ds       = 2*(prev_steering_-steering);
    prev_steering_        = steering;
    error_               += cte*cte + /*steering*steering + */ds*ds;

    const double now = time(NULL);
    if (!start_time_) {
        start_time_ = now;
        return;
    }
    
    const double diff = now-start_time_;
    if (diff < 70) { // collect error for 4 laps
        return;
    }

    // just monitoring: no learning
    std::cout << "Error: " << error_ << endl;
    start_time_ = 0;
    error_ = 0;
    return;

    // enough time (four laps) passed: update the coefficient
    history_.clear();
    prev_cte_ = 0;
    if (!is_error_set_) {
        is_error_set_ = true;
        rate_         = K/20;
        std::cout << "Initial error:" << error_ << " K: " << K << " Will try: " << K+rate_ << std::endl;
        best_error_   = error_;
        error_        = 0;
        start_time_   = now;
        K            += rate_;
        return;
    }
    
    if (is_up_) {
        if (error_ < best_error_) {
            std::cout << "New best error1:" << error_ << " K: " << K << " rate: " << rate_ << std::endl;
            best_error_ = error_;
            error_      = 0;
            start_time_ = now;
            rate_      *= 0.8;//1.1;
            K          += rate_;
            std::cout << " Will try K1: " << K << " rate " << rate_ << std::endl;
            return;
        }

        std::cout << "Error increased2:" << error_ << " K: " << K << " rate: " << rate_ << std::endl;
        is_up_      = false;
        K          -= 2*rate_;
        rate_      *= 1.1;//0.8;
        error_      = 0;
        start_time_ = now;
        std::cout << " Will try K2: " << K << " rate " << rate_ << std::endl;
        return;
    }

    is_up_ = true;
    if (error_ < best_error_) {
        std::cout << "New best error3:" << error_ << " K: " << K << " rate: " << rate_ << std::endl;
        best_error_ = error_;
        error_      = 0;
        start_time_ = now;
        rate_      *= 0.8;//1.1;
        K          += rate_;
        std::cout << " Will try K3: " << K << " rate " << rate_ << std::endl;
        return;
    }
    
    std::cout << "Error increased4:" << error_ << " K: " << K << " rate: " << rate_ << std::endl;
    K          += 2*rate_;
    rate_      *= 1.1;//0.8;
    error_      = 0;
    start_time_ = now;
    std::cout << " Will try K4: " << K << " rate " << rate_ << std::endl;
    return;
}
