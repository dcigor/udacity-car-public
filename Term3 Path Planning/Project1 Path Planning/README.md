# CarND-Path-Planning-Project
This project is based on Udacity template: https://github.com/udacity/CarND-Path-Planning-Project
   
## Architecture.
The code is split into several classes:
###Point
Represents a 2D point along with operations, such as "/" and "-".
###Points
Represents a collection of 2D points.
###Road
Represents the road, including its waypoints, tranformation from Frenet to cartesian coordinates and etc.
###OtherCar
Represents a non-controllable car on the road, provided by the sensor fusion.
###CarState
The current state of my car provided by sensor data.
###StartEnd
The beginning and end of the path trajectory. For both ends: includes current position, speed and acceleration (dot and dot_dot suffixes mean speed and acceleration correspondigly).
###Driver
The class responsible for the driver of the autonomous car.

## Implementation
###Parsing sensor data
The helper class Road parses the waypoints on startup of the program.
The helper class CarState parses JSON representation of the CarState and other sensor data at each iteration.

###High level design
The car begins driving straight in the middle line.
The car moves to the right lane whenever it is available.
When the lane is moving slowly, the car switches to the faster lane on the left to pass, then returns to the right side of the road.

###Driving in lane: Driver::keepLineXY()
We generate 1.5 seconds of waypoints for the car. The end point of the trajectory is the farthest point that can be reached given the current speed and acceleration as calculated in StartEnd::setMaxSpeed(). The starting speed and acceleration of the trajectory match the speed and trajectory of the previously generated path. The end speed matches the target speed of the road (currently 46 MPH). The end acceleration is zero.
If there is another car close to us in front (see OtherCar other = carState.nearestFrontCar()), we reduce the intended speed to match the speed of the other car to prevent the crash.
If the trajectory violates the speed and acceleration constraints, we reject the trajectory and regenerate the new one allowing more time to reach the specified end point.

###Lane changing: Driver::ChangeLane()
If there is a car close by in the target lane, we don't change the lane at this time.
Otherwise we try to generate a trajectory to change the lane in 2 seconds. If the trajectory does not satisfy the speed constraints, we try another trajectory allowing each time more time to reach the end point.
The trajectory end points are similar to above.
If it is not possible to find acceptable trajectory in 10 iterations, we don't change the lane at this time.
Once the lane change succeeds, we prevent the next line change for 10 seconds in order to allow forward progress in the new lane.

###Trajectory generation: Driver::buildTrajectory()
The trajectory is generated in Cartesian coordinates using the JMT() function from the CarND course. It generates 5-th degree polynomial that passes through the given start and end point, their speeds and accelerations and minimizes jerk (third derivative).
Then these polynomial coefficients are used to find each waypoint on the trajectory using Driver::poly() function.
