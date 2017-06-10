#include <math.h>
#include <uWS/uWS.h>

#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

#define USE_LATEST_UWS_VERSION

#ifdef USE_LATEST_UWS_VERSION
// for the latest version of UWS
#define UWS_PARAMTYPE *
#define UWS_MEMBER ->
#else
//for the version 0.13
#define UWS_PARAMTYPE
#define UWS_MEMBER .
#endif

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(const Eigen::VectorXd &xvals, const Eigen::VectorXd &yvals, const int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}
Eigen::VectorXd polyfit(const vector<double> &xvals, const vector<double> &yvals, const int order) {
    Eigen::VectorXd x(xvals.size());
    Eigen::VectorXd y(yvals.size());
    for (size_t i=0; i<xvals.size(); ++i) {
        x(i) = xvals[i];
        y(i) = yvals[i];
    }
    return polyfit(x, y, order);
}

const double Lf = 2.7;

Eigen::VectorXd globalKinematic(const double x, const double y, const double psi, const double v,
                                const double delta, const double a, const double dt) {
    Eigen::VectorXd next_state(4);
    
    // Recall the equations for the model:
    // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
    // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
    // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
    // v_[t+1] = v[t] + a[t] * dt
    next_state(0) = x   + v * cos(psi) * dt;
    next_state(1) = y   + v * sin(psi) * dt;
    next_state(2) = psi + v / Lf * delta * dt;
    next_state(3) = v   + a * dt;
    return next_state;
}

// TODO: Implement the global kinematic model.
// Return the next state.
//
// NOTE: state is [x, y, psi, v]
// NOTE: actuators is [delta, a]
Eigen::VectorXd globalKinematic(const Eigen::VectorXd &state,
                                const Eigen::VectorXd &actuators, const double dt) {
    return globalKinematic(state(0), state(1), state(2), state(3), actuators(0), actuators(1), dt);
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> UWS_PARAMTYPE ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          const vector<double> ptsx = j[1]["ptsx"];
          const vector<double> ptsy = j[1]["ptsy"];
          const double px  = j[1]["x"];
          const double py  = j[1]["y"];
          const double psi = j[1]["psi"];
          const double v   = j[1]["speed"];

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */

            //predict where the car is going to be in 100 ms in global coordinates
            const Eigen::VectorXd globalPredictedPos = globalKinematic(px, py, psi, v, 0, 0, 0);
            const double globalPredictedX = globalPredictedPos(0);
            const double globalPredictedY = globalPredictedPos(1);
            
            const World w(globalPredictedX, globalPredictedY, psi);
            const vector<Eigen::Vector3d> pts_local(w.globalToLocal(ptsx, ptsy));
            //const Eigen::Vector3d localPredictPos    = w.globalToLocal(globalPredictedPos(0), globalPredictedPos(1));
            
          json msgJson;

            //Display the MPC predicted trajectory
            //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
            // the points in the simulator are connected by a Green line
            const vector<double> mpc_x_vals(w.getRow(0, pts_local));
            const vector<double> mpc_y_vals(w.getRow(1, pts_local));
            msgJson["mpc_x"] = mpc_x_vals;
            msgJson["mpc_y"] = mpc_y_vals;
            
            // fit the 3rd degree polynomial, in vehicle coordinate system
            const Eigen::VectorXd pol = polyfit(mpc_x_vals, mpc_y_vals, 3);

            // calculate cte as distance from the current vehicle position to the polynomial
            const double polynomial_y = polyeval(pol, 0);

            const double cte  = polynomial_y;
            const double epsi = -atan(pol[1]);
            cout << "POLY " << epsi << " " << endl;
            
            mpc.Solve(0, 0, 0, v, cte, epsi, pol);
            const double throttle_value = mpc.getAcceleration();
            const double steer_value    = -mpc.getSteeringDelta();
            
            // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
            // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
            msgJson["steering_angle"] = steer_value;
            msgJson["throttle"      ] = throttle_value;

          //Display the waypoints/reference line
            //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
            // the points in the simulator are connected by a Yellow line
            vector<double> next_x_vals;
            vector<double> next_y_vals;
            
            // take into account the steering
            const World w2(globalPredictedX, globalPredictedY, psi+steer_value/Lf/2);
            const vector<Eigen::Vector3d> pts_local2(w2.globalToLocal(ptsx, ptsy));
            next_x_vals = w2.getRow(0, pts_local2);
            next_y_vals = w2.getRow(1, pts_local2);

            msgJson["next_x"] = mpc.getRow(0);
            msgJson["next_y"] = mpc.getRow(1);

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          //this_thread::sleep_for(chrono::milliseconds(100));
          ws UWS_MEMBER send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws UWS_MEMBER send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> UWS_PARAMTYPE ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> UWS_PARAMTYPE ws, int code,
                         char *message, size_t length) {
    ws UWS_MEMBER close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
