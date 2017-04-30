#ifndef PID_H
#define PID_H

#include <vector>

class PID {
public:
  /*
  * Errors
  */
  double p_error;
  double i_error;
  double d_error;

  /*
  * Coefficients
  */ 
  double Kp;
  double Ki;
  double Kd;

  /*
  * Constructor
  */
  PID();

  /*
  * Destructor.
  */
  virtual ~PID();

  /*
  * Initialize PID.
  */
  void Init(double Kp, double Ki, double Kd);

  /*
  * Update the PID error variables given cross track error.
  */
  void UpdateError(double cte);

  /*
  * Calculate the total PID error.
  */
  double TotalError();
private:
    double prev_cte_ = 0;
    std::vector<double> history_;

    // learn coefficients
    void learn_coefficient(double cte, double &K);
    double best_error_ = 0, error_ = 0;
    bool is_error_set_ = false;
    int start_time_    = 0;
    double rate_       = 1;
    bool is_up_        = true;
    double prev_steering_ = 0;
};

#endif /* PID_H */
