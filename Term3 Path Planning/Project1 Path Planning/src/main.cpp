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

    // remove the elements from the end to reach the desired size
    void trim(const size_t size) {
        while (x.size()>size) {
            x.pop_back();
            y.pop_back();
        }
    }
    
    friend ostream& operator<<(ostream& os, const Points& points) {
        os << "[";
        for (size_t i=0; i<points.size(); ++i) {
            os << points[i].x;
        }
        os << "]";
        return os;
    }
};

class Road {
public:
    static double mph2mps(const double mph) {
        return mph * 1609 / 3600;
    }
    static double mps2mph(const double mps) {
        return mps / 1609 * 3600;
    }
    const double TARGET_SPEED_MPH = 45;
    const double TARGET_SPEED     = mph2mps(TARGET_SPEED_MPH);

    const double LANE_WIDTH = 4;

    double lineCenter(const int lane) const {
        return LANE_WIDTH*lane - LANE_WIDTH/2;
    }

    // loads waypoints on start
    Road(const string &filename, const double in_max_s) : max_s(in_max_s) {
        // Waypoint map to read from
        ifstream in_map_(filename.c_str(), ifstream::in);

        string line;
        while (getline(in_map_, line)) {
            istringstream iss(line);
            double x, y;
            float s, d_x, d_y;
            iss >> x >> y >> s >> d_x >> d_y;
            map_waypoints_x_ .push_back(x);
            map_waypoints_y_ .push_back(y);
            map_waypoints_s_ .push_back(s);
            map_waypoints_dx_.push_back(d_x);
            map_waypoints_dy_.push_back(d_y);
        }
    }

    Point getXY(const double s, const double d) const {
        return Point::fromArray(::getXY(s, d, map_waypoints_s_, map_waypoints_x_, map_waypoints_y_));
    }

    const double max_s;
private:
    vector<double> map_waypoints_x_;
    vector<double> map_waypoints_y_;
    vector<double> map_waypoints_s_;
    vector<double> map_waypoints_dx_;
    vector<double> map_waypoints_dy_;
};

struct OtherCar {
    int id;
    Point xy;
    Point speedXY;
    double s, d;
    
    static OtherCar initFromVector(const vector<double> &v) {
        OtherCar car;
        car.id        = v[0];
        car.xy.x      = v[1];
        car.xy.y      = v[2];
        car.speedXY.x = v[3];
        car.speedXY.y = v[4];
        car.s         = v[5];
        car.d         = v[6];
        return car;
    }
    
    static OtherCar createInvalid() {
        OtherCar car;
        car.id = -1;
        return car;
    }
    
    bool isValid() const {return id >= 0;}
    
    double predictS(const double t) const {
        // assume the car goes straight in the lane
        return s + t*speedXY.abs();
    }
};

struct CarState {
    Point  xy;
    double s;
    double d;
    double yaw;
    double speed;
    
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
        carState.xy.x  = j[1]["x"];
        carState.xy.y  = j[1]["y"];
        carState.s     = j[1]["s"];
        carState.d     = j[1]["d"];
        carState.yaw   = double(j[1]["yaw"])/360*2*M_PI;
        carState.speed = Road::mph2mps(j[1]["speed"]);

        const vector<double> previous_path_x = j[1]["previous_path_x"];
        const vector<double> previous_path_y = j[1]["previous_path_y"];
        carState.previous_pathXY.x           = previous_path_x;
        carState.previous_pathXY.y           = previous_path_y;
        carState.end_path_s                  = j[1]["end_path_s"];
        carState.end_path_d                  = j[1]["end_path_d"];

        const vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];
        carState.sensor_fusion = sensor_fusion;
        return carState;
    }
    
    OtherCar nearestFrontCar(const int lane, const double s, const Road &road) const {
        OtherCar other = OtherCar::createInvalid();
        for (const auto v : sensor_fusion) {
            const OtherCar curr = OtherCar::initFromVector(v);
            if (curr.s < s) {
                continue; // this car is behind. ignore it for now
            }
            
            const double lineCenter = road.lineCenter(lane);
            if (abs(curr.d-lineCenter) > road.LANE_WIDTH/2) {
                continue; // this car is in another line - ignore it
            }
            
            if (!other.isValid() || curr.s < other.s) {
                other = curr;
            }
        }
        return other;
    }

    friend ostream& operator<<(ostream& os, const CarState& s) {
        os << "{" << s.xy << ", s: " << s.s << ", d: " << s.d << ", ya: " << s.yaw << ", sp: " << s.speed << ", end s: " <<
            s.end_path_s << ", d: " << s.end_path_d << "}";
        return os;
    }
};

class StartEnd {
public:
    StartEnd(const Road &road, const CarState &carState, const Points &points, const double in_TICK, const double in_endD) : TICK_(in_TICK) {
        startS          = carState.s;
        startD          = carState.d;
        startXY         = carState.xy;
        startXY_dot     = {carState.speed * cos(carState.yaw), carState.speed * sin(carState.yaw)};
        startXY_dot_dot = {0, 0};
        desiredEndSpeed = 1000000;

        if (points.size()>2) {
            startS                = carState.end_path_s;
            startD                = carState.end_path_d;
            startXY               = points[points.size()-1];
            const Point prev      = points[points.size()-2];
            startXY_dot           = (startXY-prev)/TICK_;

            const Point prev_prev = points[points.size()-3];
            const Point prevSpeed = (prev-prev_prev)/TICK_;
            startXY_dot_dot       = (startXY_dot-prevSpeed)/TICK_;
            
            if ((road.getXY(startS, startD)-startXY).abs() > 0.01) {
                cout << "ERR " << (road.getXY(startS, startD)-startXY).abs() << endl;
            }
        }

        endD = in_endD;
    }

    void setMaxSpeed(const Road &road, const double maxSpeed, const double maxAcc, const double T) {
        desiredEndSpeed = min(desiredEndSpeed, maxSpeed); // keep reducing speed with each additional restriction
        endS = startS + desiredEndSpeed * T;

        // can we reach the destination within the max acceleration?
        const double maxS = startS + startXY_dot.abs() * T + 0.5*maxAcc * T * T;
        if (endS > maxS) {
            endS = maxS;
            // reduce destination speed
            desiredEndSpeed = min(maxSpeed, startXY_dot.abs()+maxAcc * T);
        }
        const double minS = startS + startXY_dot.abs() * T - 0.5*maxAcc * T * T;
        if (endS < minS) {
            endS = minS;
            // increase destination speed, as cannot slow down that fast
            desiredEndSpeed = max(desiredEndSpeed, startXY_dot.abs()-maxAcc * T);
        }

        setMaxEndS(road, endS);
    }

    void setMaxEndS(const Road &road, const double in_endS) {
        if (in_endS <= endS) {
            endS = in_endS;
            const double endS_prev = endS - desiredEndSpeed * TICK_; // the point before the last point in the future

            endXY = road.getXY(endS, endD);
            const auto endXY_prev = road.getXY(endS_prev, endD);
            const double endYaw   = atan2(endXY.y-endXY_prev.y, endXY.x-endXY_prev.x);
            
            endXY_dot  = {desiredEndSpeed * cos(endYaw), desiredEndSpeed * sin(endYaw)};
        }
    }

    Point startXY;
    Point startXY_dot;
    Point startXY_dot_dot = {0, 0};
    Point endXY, endXY_dot, endXY_dot_dot = {0,0};

    double startS, startD;
    double desiredEndSpeed;
    double endS, endD;
private:
    const double TICK_;
};

class Driver {
public:
    const double TICK = 0.02; // 20 ms between points

    Driver(const Road &road) : road_(road) {}

    Points drive(const CarState& carState) {
        const auto now = std::chrono::steady_clock::now();
        if ((carState.speed>10) && (std::chrono::duration_cast<std::chrono::seconds>(now-safe_time_).count() > 0)) {
            if (lane_ != 3) {
                const Points points = changeLine(carState, lane_+1);
                if (points.size() > 0) {
                    return points;
                }
            }

            if ((lane_ != 1) && (carState.speed < road_.TARGET_SPEED-3)) {
                const Points points = changeLine(carState, lane_-1);
                if (points.size() > 0) {
                    return points;
                }
            }
        }
        return keepLineXY(carState, lane_);
    }

    Points keepLineXY(const CarState& carState, const int lane) {
        // go from the current state to state with the same d, s+10, speed : 10
        const double T       = 1.5;    // reach the target speed within this time
        const int count      = T/TICK; // how many points?
        const double MAX_ACC = 3;      // do not exceed this acceleration

        // use previously generated points, when present
        Points points = carState.previous_pathXY;
        if (points.size()>count/2) {
            ///cout << "("<<points.size()<<") ";
            return points; // keep following the earlier generated path
        }

        StartEnd startEnd(road_, carState, points, TICK, road_.lineCenter(lane));

        // can we reach the destination within the max acceleration?
        startEnd.setMaxSpeed(road_, road_.TARGET_SPEED, MAX_ACC, T);

        // reduce the speed if there is a car in front
        const OtherCar other = carState.nearestFrontCar(lane, carState.s, road_);
        if (other.isValid()) {
            const double SAFE_TIME = 4;
            const double s = other.predictS(T-SAFE_TIME);
            if (s < startEnd.endS) {
                startEnd.setMaxSpeed(road_, (startEnd.startXY_dot.abs()+other.speedXY.abs())/2-2, MAX_ACC, T);
            }
/*                cout << "[D " << startEnd.endS - s << "] -> " << Road::mps2mph((other.speedXY.abs()-2)) << endl;
                // there is another car within N seconds, go a little slower than it
                startEnd.setMaxSpeed(road_, other.speedXY.abs()-2, MAX_ACC, T);
            }   else {
                const double s = other.predictS(T-SAFE_TIME-1);
                if (s < startEnd.endS) {
                    cout << "[~ " << startEnd.endS - s << "] -> " << Road::mps2mph((startEnd.startXY_dot.abs()+other.speedXY.abs())/2) << endl;
                    // there is another car within N+1 seconds, match its speed
                    startEnd.setMaxSpeed(road_, (startEnd.startXY_dot.abs()+other.speedXY.abs())/2-2, MAX_ACC, T);
                }
            }
  */      }

        return buildTrajectory(points, startEnd, T);
    }

    Points changeLine(const CarState& carState, const int lane) {
        const double T       = 3; // reach the target speed within this time
        const double MAX_ACC = 1; // do not exceed this acceleration

        Points points = carState.previous_pathXY;
        // shorten previous path
        //points.trim(20);

        StartEnd startEnd(road_, carState, points, TICK, road_.lineCenter(lane));
        //startEnd.setMaxSpeed(road_, startEnd.startXY_dot.abs(), MAX_ACC, T);
        startEnd.setMaxSpeed(road_, road_.TARGET_SPEED, MAX_ACC, T);
        startEnd.setMaxEndS (road_, startEnd.endS-road_.LANE_WIDTH/2);

        const OtherCar other = carState.nearestFrontCar(lane, carState.s-20, road_);
        if (other.isValid()) {
            const double SAFE_TIME = 3;
            const double s = other.predictS(T-SAFE_TIME);
            if (s < startEnd.endS) {
                cout << "cant change to " << lane << endl;
                return {}; // another car is too close
            }
        }

        lane_      = lane;
        safe_time_ = std::chrono::steady_clock::now() + std::chrono::seconds(10);
        return buildTrajectory(points, startEnd, T);
    }

    Points buildTrajectory(Points points, const StartEnd &startEnd, const double T) {
        const int count    = T/TICK;
        const auto coefs_x = JMT({startEnd.startXY.x, startEnd.startXY_dot.x, startEnd.startXY_dot_dot.x}, {startEnd.endXY.x, startEnd.endXY_dot.x, startEnd.endXY_dot_dot.x}, T);
        const auto coefs_y = JMT({startEnd.startXY.y, startEnd.startXY_dot.y, startEnd.startXY_dot_dot.y}, {startEnd.endXY.y, startEnd.endXY_dot.y, startEnd.endXY_dot_dot.y}, T);
            
        string s = to_string(coefs_x[0]) + "+" +        to_string(coefs_x[1]) + "*x + "      + to_string(coefs_x[2]) + "*x*x + " +
                       to_string(coefs_x[3]) + "*x*x*x + "+ to_string(coefs_x[4]) + "*x*x*x*x +" + to_string(coefs_x[5]) + "*x*x*x*x*x";

        cout << s << endl;
        cout << "START s:" << startEnd.startS << ", " << startEnd.startXY << ", Speed: " << startEnd.startXY_dot << ", A: " << startEnd.startXY_dot_dot << ", " << startEnd.startXY_dot_dot.abs() << ", end: " << startEnd.endXY << endl;
        cout << "D S: " << startEnd.endS - startEnd.startS << ", dxy: " << (startEnd.endXY-startEnd.startXY).abs() << ", speed: " << (startEnd.endS - startEnd.startS)/T << ", mph: " << Road::mps2mph((startEnd.endS - startEnd.startS)/T) << endl;
        cout << "PT: " << road_.getXY(startEnd.startS, startEnd.startD) << " -> " << road_.getXY(startEnd.endS, startEnd.endD) << ", D: " << (road_.getXY(startEnd.startS, startEnd.startD)-road_.getXY(startEnd.endS, startEnd.endD)).abs() << endl;

        Point last      = startEnd.startXY;
        Point lastSpeed = startEnd.startXY_dot;
        for (double t   = TICK; t < T && points.size()<count*2; ) {
            const Point curr = {poly(coefs_x, t), poly(coefs_y, t) };
            const Point speed = (curr-last) / TICK;
            const Point acc   = (speed-lastSpeed) / TICK;
            if ((speed.abs() > road_.TARGET_SPEED+0.4)/* ||
                (acc.abs() > 9)*/) {
                cout << "TOO FAST: " << curr << ", Speed: " << speed << " mps, " << road_.mps2mph(speed.abs()) << " mph, acc: " << acc << endl;
                //t -= TICK/200; // reduce the step, so we can reach it with smaller speed
                //continue;
            }

            points.push(curr);
            cout << "POINT: " << curr << ", Speed: " << speed << " mps, " << road_.mps2mph(speed.abs()) << " mph, acc: " << acc << endl;
            last      = curr;
            lastSpeed = speed;
            t        += TICK;
        }
        return points;
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
    Road road_;
    int lane_ = 2;
    std::chrono::time_point<std::chrono::steady_clock> safe_time_ = std::chrono::steady_clock::now();
};

int main() {
    uWS::Hub h;

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    Road road("../data/highway_map.csv", 6945.554);
    Driver driver(road);

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
            const CarState carState = CarState::initFromJson(j);
            const Points points     = driver.drive(carState);

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

