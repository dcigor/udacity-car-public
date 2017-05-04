#include "PID.h"
#include <time.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>

/*
 Reflections: 
 
 The P coefficient is proportional to the Cross Track Error (CTE).
 Larger values of P return the car to the middle of the track faster. 
 This works OK for driving on a straight segment of the track.
 However, turning segment of the track might present a problem. 
 Specifically: with large P value, the car crosses the center of the track at a 
 significant angle. If the track turns at this point into the opposite side, 
 the car might be facing the end of the track at even steeper angle.
 It might be impossible to steer back to the center of the road from this position.
 
 From another hand, small values of P do not allow the car to follow the curve 
 of the track when the track turns sharply.
 
 Using the learn_coefficient() method, it is possible to find the most stable range for 
 the P coefficient. For throttle 0.4, the good range for P is 0.1-0.3. 
 By continuing using learn_coefficient() method, the best P cefficient comes to 0.25.
 
 ------
 
 The D coefficient reduces oversteering. With small values of D, the car oscillates around 
 the center of the track. It causes two problems: first unpleasant experience. Second:
 at a sharp turn of the track, it increases the chance that the car will be facing the opposite 
 end of the track during one of the oscillations. In some of these occasions, the angle might be too 
 steep for the car to recover.
 
 Large values of D do not allow the car enough freedom to recover when it finds itself in a situation
 when a sharp turn is required. Large D will force the car to go straight instead.
 
 Again, learn_coefficient() method suggests the best value for the D coefficient is about 8.
 
 ------
 
 The I coefficient allows correction of the bias. For example, if the car believes it goes straight,
 when in fact its wheels point to a side.
 
 It appears that the simulator does not have the bias. So, although I added a way to use 
 the I coefficient by integrating CTE over hundred(s) of samples,
 for solving the bias, it is not yet possible to test if this implementation is the best one.
 
 --------
 
 Future improvements:
 It appears that the PID algorithm works best when the track is straight, because it 
 tries to put the car into the middle of the track looking straight.
 
 However, there are two problems with that.
 1. When the track has sharp turn, the best car trajectory is not the middle of the track.
 Instead the car better follow a curve with the biggest radius (within the allowed portion of the road) as the racecars do on TV.
 
 2. More importantly, the best heading of the car during a turn is not pointing forward.
 Instead, the best heading is turning at an angle that follows the turning angle.
 However, PID algorithm instead will try to point the car forward.
 
 In order to solve the second issue, we would have find the turning angle of the track in front of the car
 using for example Line Detection, then center the PID around this angle. 
 
 Basically, in calculation of d_error, instead of 
    d_error = cte-prev_cte_;
 Which tries to point the car straight, we need to calculate d_error by pointing the car somewhat invards.
 
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
    if (history_.size() > 200) {
        history_.erase(history_.begin());
    }

    p_error = cte;
    d_error = cte-prev_cte_;
    i_error = std::accumulate(history_.begin(), history_.end(), 0);

    prev_cte_ = cte;

    // this allow to learn the best value of any of the coefficients Kp, Kd or Ki
    //learn_coefficient(cte, Kp);
}

double PID::TotalError() {
    const double res = -(p_error * Kp + d_error * Kd + i_error * Ki);
    return std::min(std::max(res, -1.), 1.);
}

// this is only needed to find the best values of the coefficients.
void PID::learn_coefficient(const double cte, double &K) {
    const double steering = TotalError();
    const double ds       = 10000*(prev_steering_-steering);
    prev_steering_        = steering;
    error_               += cte*cte*10000 + /*steering*steering + */ds*ds; // choose what needs to be minimized

    const double now = time(NULL);
    if (!start_time_) {
        start_time_ = now;
        return;
    }
    
    const double diff = now-start_time_;
    if (diff < 70*4) { // collect error for 4 laps
        return;
    }

    // enough time (four laps) passed: update the coefficient
    history_.clear();
    prev_cte_ = 0;
    if (!is_error_set_) {
        is_error_set_ = true;
        rate_         = K/20;
        std::cout << "Initial error: " << error_ << " K: " << K << " Will try: " << K+rate_ << std::endl;
        best_error_   = error_;
        error_        = 0;
        start_time_   = now;
        K            += rate_;
        return;
    }
    
    if (is_up_) {
        if (error_ < best_error_) {
            std::cout << "New best error1: " << error_ << " K: " << K << " rate: " << rate_ << std::endl;
            best_error_ = error_;
            error_      = 0;
            start_time_ = now;
            rate_      *= 0.8;
            K          += rate_;
            std::cout << " Will try K1: " << K << " rate " << rate_ << std::endl;
            return;
        }

        std::cout << "Error increased2: " << error_ << " K: " << K << " rate: " << rate_ << std::endl;
        is_up_      = false;
        K          -= 2*rate_;
        rate_      *= 1.1;
        error_      = 0;
        start_time_ = now;
        std::cout << " Will try K2: " << K << " rate " << rate_ << std::endl;
        return;
    }

    is_up_ = true;
    if (error_ < best_error_) {
        std::cout << "New best error3: " << error_ << " K: " << K << " rate: " << rate_ << std::endl;
        best_error_ = error_;
        error_      = 0;
        start_time_ = now;
        rate_      *= 0.8;
        K          += rate_;
        std::cout << " Will try K3: " << K << " rate " << rate_ << std::endl;
        return;
    }

    std::cout << "Error increased4: " << error_ << " K: " << K << " rate: " << rate_ << std::endl;
    K          += 2*rate_;
    rate_      *= 1.1;
    error_      = 0;
    start_time_ = now;
    std::cout << " Will try K4: " << K << " rate " << rate_ << std::endl;
}
