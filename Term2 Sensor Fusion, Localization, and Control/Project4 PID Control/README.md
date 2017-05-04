# CarND-Controls-PID
This project is based on Udacity template: https://github.com/udacity/CarND-PID-Control-Project

##Reflections: 
 
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
 
##Future improvements:
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