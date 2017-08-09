#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
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

using namespace std;
using namespace Eigen;

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
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(const double s, const double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

    while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) )) {
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

vector<double> JMT(const vector<double> &start, const vector<double> &end, const double T) {
    /*
     Calculate the Jerk Minimizing Trajectory that connects the initial state
     to the final state in time T.
     
     INPUTS
     
     start - the vehicles start location given as a length three array
     corresponding to initial values of [s, s_dot, s_double_dot]
     
     end   - the desired end state for vehicle. Like "start" this is a
     length three array.
     
     T     - The duration, in seconds, over which this maneuver should occur.
     
     OUTPUT
     an array of length 6, each value corresponding to a coefficent in the polynomial
     s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5
     
     EXAMPLE
     
     > JMT( [0, 10, 0], [10, 10, 0], 1)
     [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
     */
    
    const double si = start[0], si_d = start[1], si_dd = start[2];
    const double sf = end  [0], sf_d =   end[1], sf_dd =   end[2];
    Matrix3f A;
    Vector3f b;
    A << T*T*T, T*T*T*T, T*T*T*T*T,
    3*T*T, 4*T*T*T, 5*T*T*T*T,
    6*T  , 12 *T*T, 20 *T*T*T;
    b[0] = sf    - (si+si_d*T+0.5*si_dd*T*T);
    b[1] = sf_d  - (si_d + si_dd*T);
    b[2] = sf_dd - si_dd;

    const Vector3f x = A.colPivHouseholderQr().solve(b);
    return {si, si_d, 0.5*si_dd, x[0], x[1], x[2]};
}

class Driver {
public:
    static double mph2mps(const double mph) {
        return mph * 1609 / 3600;
    }
    static double mps2mph(const double mps) {
        return mps / 1609 * 3600;
    }
    const double TARGET_SPEED_MPH = 45;
    const double TARGET_SPEED = mph2mps(TARGET_SPEED_MPH);
    const double TICK = 0.02; // 20 ms between points

    struct Point {
        double x, y;

        static Point fromArray(const vector<double> &a) {
            Point p;
            p.x = a[0];
            p.y = a[1];
            return p;
        }
        
        Point operator-(const Point &b) const {
            return {x-b.x, y-b.y};
        }
        
        Point operator/(const double &b) const {
            return {x/b, y/b};
        }

        double abs() const {
            return sqrt(x*x+y*y);
        }
        
        friend ostream& operator<<(ostream& os, const Point& p) {
            os << "{" << p.x << ", " << p.y << ", |" << p.abs() <<"|}";
            return os;
        }
    };

    struct Points {
        vector<double> x;
        vector<double> y;

        void push(const Point &point) {
            x.push_back(point.x);
            y.push_back(point.y);
        }
        
        size_t size() const {
            return x.size();
        }

        Point operator[](const size_t i) const {
            return {x[i], y[i]};
        }
        
        Points& operator+=(const Points& b) {
            x.insert(x.end(), b.x.begin(), b.x.end());
            y.insert(y.end(), b.y.begin(), b.y.end());
            return *this;
        }
    };

    struct CarState {
        Point  carXY;
        double car_s;
        double car_d;
        double car_yaw;
        double car_speed;

        // Previous path data given to the Planner
        Points previous_pathXY;
        // Previous path's end s and d values
        double end_path_s;
        double end_path_d;
        
        // Sensor Fusion Data, a list of all other cars on the same side of the road.
        vector<vector<double>> sensor_fusion;

        template<typename J>
        static CarState initFromJson(const J& j) {
            CarState carState;
            carState.carXY.x   = j[1]["x"];
            carState.carXY.y   = j[1]["y"];
            carState.car_s     = j[1]["s"];
            carState.car_d     = j[1]["d"];
            carState.car_yaw   = double(j[1]["yaw"])/360*2*M_PI;
            carState.car_speed = j[1]["speed"];

            const vector<double> previous_path_x = j[1]["previous_path_x"];
            const vector<double> previous_path_y = j[1]["previous_path_y"];
            carState.previous_pathXY.x = previous_path_x;
            carState.previous_pathXY.y = previous_path_y;
            carState.end_path_s        = j[1]["end_path_s"];
            carState.end_path_d        = j[1]["end_path_d"];

            const vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];
            carState.sensor_fusion = sensor_fusion;
            return carState;
        }

        friend ostream& operator<<(ostream& os, const CarState& s) {
            os << s.carXY << ", s: " << s.car_s << ", d: " << s.car_d << ", ya: " << s.car_yaw << ", sp: " << s.car_speed;
            return os;
        }
    };

    // loads waypoints on start
    Driver(const string &filename) {
        // Waypoint map to read from
        ifstream in_map_(filename.c_str(), ifstream::in);

        string line;
        while (getline(in_map_, line)) {
            istringstream iss(line);
            double x, y;
            float s, d_x, d_y;
            iss >> x >> y >> s >> d_x >> d_y;
            map_waypoints_x .push_back(x);
            map_waypoints_y .push_back(y);
            map_waypoints_s .push_back(s);
            map_waypoints_dx.push_back(d_x);
            map_waypoints_dy.push_back(d_y);
        }
    }

    Points drive(const CarState& carState) {
        return keepLineXY(carState);
    }

    Points keepLine(const CarState& carState, const int line=2) {
        return keepLineXY(carState);
    
        // go from the current state to state with the same d, s+10, speed : 10
        const double T    = 1; // reach the target speed within this time
        const int count   = 50; // how many points?
        const double tick = T/count;

        // end state in Frenet
        const double startS = carState.car_s;
        const double startS_dot = carState.car_speed;
        const double startS_dot_dot = 0;
        const double endS = carState.car_s + TARGET_SPEED * T;
        const double endS_dot = TARGET_SPEED;
        const double endS_dot_dot = 0;

        const double startD = carState.car_d;
        const double endD = carState.car_d;
        const double endD_dot = 0;
        const double endD_dot_dot = 0;

        const auto coefs_s = JMT({startS, startS_dot, startS_dot_dot}, {endS, endS_dot, endS_dot_dot}, T);
        const auto coefs_d = JMT({startD, 0, 0}, {endD, 0, 0}, T);
    
        cout << "AA ";
        Points points;
        double lastS = startS;
        for (int i=0; i<count; ++i) {
            const double s = poly(coefs_s, i*tick);
            const double d = startD;//poly(coefs_d, T/count*i);
            if (s<startS || s>endS) {
                cout << s;
            }
            cout << endS-s << ", speed: " << (s-lastS)/tick << endl;
            lastS = s;

            const Point xy = getXY(s, d);
            points.x.push_back(xy.x);
            points.y.push_back(xy.y);
        }
        cout << "BB ";
        return points;
    }

    Points keepLineXY(const CarState& carState, const int lane=2) {
        // go from the current state to state with the same d, s+10, speed : 10
        const double T    = 2;      // reach the target speed within this time
        const int count   = T/TICK; // how many points?
        const double MAX_ACC = 5;   // do not exceed this acceleration

        Point startXY     = carState.carXY;
        Point startXY_dot = {carState.car_speed * cos(carState.car_yaw),
                             carState.car_speed * sin(carState.car_yaw)};
        Point startXY_dot_dot = {0, 0};

        // use previously generated points, when present
        Points points = carState.previous_pathXY;
        if (points.size()>count/2) {
            return points; // keep following the earlier generated path
        }

        double startS = carState.car_s;
        if (points.size()) {
            startXY     = points[points.size()-1];
            const Point prev = points[points.size()-2];
            startXY_dot = (startXY-prev)/TICK;
            const Point prev_prev = points[points.size()-3];
            const Point prevSpeed = (prev-prev_prev)/TICK;
            startXY_dot_dot = (startXY_dot-prevSpeed)/TICK;

            startS = carState.end_path_s;
        }

        double desiredEndSpeed = TARGET_SPEED;
        double endS = startS + TARGET_SPEED * T, endD = lane * 4 - 2;

        // how far we can go with the max acceleration?
        const double maxS = startS + startXY_dot.abs() * T + 0.5*MAX_ACC * T * T;
        const bool useConstantAcceleration = endS > maxS;
        if (useConstantAcceleration) {
            endS = maxS;
            desiredEndSpeed = min(TARGET_SPEED, startXY_dot.abs()+MAX_ACC * T);
        }

        const double endS_prev = endS - desiredEndSpeed * TICK; // the point before the last point in the future

        const auto endXY       = getXY(endS     , endD);
        const auto endXY_prev  = getXY(endS_prev, endD);

        const double endYaw    = atan2(endXY.y-endXY_prev.y, endXY.x-endXY_prev.x);

        const Point endXY_dot  = {desiredEndSpeed * cos(endYaw), desiredEndSpeed * sin(endYaw)};

        const Point endXY_dot_dot = {0, 0};
        
        cout << carState << ", A: " << startXY_dot_dot << ", " << startXY_dot_dot.abs() << ", end: " << endXY << endl;

        //if (useConstantAcceleration) {
        //    return points += constantAcceleration(startXY, startXY_dot, startXY_dot_dot, endXY, endXY_dot, endXY_dot_dot, T);
        //}

        const auto coefs_x = JMT({startXY.x, startXY_dot.x, startXY_dot_dot.x}, {endXY.x, endXY_dot.x, endXY_dot_dot.x}, T);
        const auto coefs_y = JMT({startXY.y, startXY_dot.y, startXY_dot_dot.y}, {endXY.y, endXY_dot.y, endXY_dot_dot.y}, T);

        Point last = startXY;
        Point lastSpeed = startXY_dot;
        for (size_t i=1; i<=count; ++i) {
            const Point curr = {poly(coefs_x, i*TICK), poly(coefs_y, i*TICK) };
            if (points.size()<count*2) {
                points.push(curr);

                const Point speed = (curr-last) / TICK;
                const Point acc   = (speed-lastSpeed) / TICK;

                cout << "POINT: " << curr << ", Speed: " << speed << " mps, " << mps2mph(speed.abs()) << " mph, acc: " << acc << endl;
                last      = curr;
                lastSpeed = speed;
            }
        }
        return points;
    }

    Points keepLineFrenet(const CarState& carState, const int lane=2) {
        // go from the current state to state with the same d, s+10, speed : 10
        const double T    = 1;      // reach the target speed within this time
        const int count   = T/TICK; // how many points?
        
        double startS = carState.car_s;
        double startD = carState.car_d;
        
        double startS_dot = carState.car_speed;
        double startD_dot = 0;
        
        double startS_dot_dot = 0;
        double startD_dot_dot = 0;

        // use previously generated points, when present
        Points points = carState.previous_pathXY;
        if (points.size()) {
            startS     = carState.end_path_s;
            startD     = carState.end_path_d;
            startS_dot = (points[points.size()-1]-points[points.size()-2]).abs()/TICK;// assume 100% of speed goes to s
        }

        const double endS = startS + TARGET_SPEED*T, endD = lane * 4 - 2;

        const double endS_dot = TARGET_SPEED, endD_dot = 0;

        const double endS_dot_dot = 0, endD_dot_dot = 0;

        const auto coefs_s = JMT({startS, startS_dot, startS_dot_dot}, {endS, endS_dot, endS_dot_dot}, T);
        const auto coefs_d = JMT({startD, startD_dot, startD_dot_dot}, {endD, endD_dot, endD_dot_dot}, T);
        
        cout << "STATE: " << carState << ", Start: " << getXY(startS, startD) << startS << ", " << startD << ", end: " << endS << ", sp: " << startS_dot << endl;

        if (points.size()<count/2) {
            Point last      = getXY(startS, startD);
            Point lastSpeed = {startS_dot, startD_dot};
            for (size_t i=1; i<=count; ++i) {
                const double currS = poly(coefs_s, i*TICK);
                const double currD = poly(coefs_d, i*TICK);
                const Point curr   = getXY(currS, currD);
                if (points.size()<count*2) {
                    points.push(curr);

                    const Point speed = (curr-last) / TICK;
                    const Point acc   = (speed-lastSpeed) / TICK;

                    cout << "POINT: " << curr << ", Speed: " << speed << " mps, " << mps2mph(speed.abs()) << " mph, acc: " << acc << endl;
                    last      = curr;
                    lastSpeed = speed;
                }
            }
        }
        return points;
    }

    
    Points keepLineSimple(const CarState& carState, const int lane=2) { // 1: left lane; 2: middle lane; 3: right lane
        // go from the current state to state with the same d, s+10, speed : 10
        const double T      = 1; // reach the target speed within this time
        const int count     = T/TICK; // how many points?
        const double startS = carState.car_s;
        const double startD = lane * 4 - 2;

        Points points;
        for (int i=0; i<count; ++i) {
            const double s = startS + TARGET_SPEED*i*TICK;
            points.push(getXY(s, startD));
        }
        return points;
    }

    Points driveCircle(const CarState& carState) {
        Points points = carState.previous_pathXY;
        Point  pos    = carState.carXY;
        double angle  = carState.car_yaw;
        const size_t path_size = carState.previous_pathXY.size();
        if (path_size) {
            pos = carState.previous_pathXY[path_size-1];
            const Point pos2 = carState.previous_pathXY[path_size-2];
            angle = atan2(pos.y-pos2.y, pos.x-pos2.x);
        }

        double dist_inc = 0.5;
        for (size_t i=0; points.size()<25/*i < 50-path_size*/; i++) {
            pos.x += dist_inc*cos(angle+(i+1)*(pi()/100));
            pos.y += dist_inc*sin(angle+(i+1)*(pi()/100));
            points.push(pos);

            if (i) {
                const double speed = (pos-points[points.size()-2]).abs() / TICK;
                cout << pos.x << ", " << pos.y << ", Speed: " << speed << "mps, " << mps2mph(speed) << "mph" << endl;
            }
        }

        return points;
    }
    
    Points constantAcceleration(const Point &startXY, const Point &startXY_dot, const Point &startXY_dot_dot, const Point &endXY, const Point &endXY_dot, const Point &endXY_dot_dot, const double T) {
        Points points;
        for (double t = 0; t<T+TICK/2; t += TICK) {
            //points.x.push_back();
        }
        return points;
    }

    Point getXY(const double s, const double d) const {
        return Point::fromArray(::getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y));
    }

    double poly(const vector<double> &coeffs, const double x) const {
        double res = 0, arg = 1;
        for (const double coef : coeffs) {
            res += coef * arg;
            arg *= x;
        }
        return res;
    }
private:
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;
};


int main() {
    uWS::Hub h;

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    Driver driver("../data/highway_map.csv");
    const double max_s = 6945.554;

    h.onMessage([&driver](uWS::WebSocket<uWS::SERVER> UWS_PARAMTYPE ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
            const Driver::CarState carState = Driver::CarState::initFromJson(j);
            const Driver::Points points     = driver.drive(carState);

            json msgJson;
          	msgJson["next_x"] = points.x;
          	msgJson["next_y"] = points.y;

          	const auto msg = "42[\"control\","+ msgJson.dump()+"]";
            //cout << msg << endl;

          	//this_thread::sleep_for(chrono::milliseconds(1000));
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
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
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
















































































