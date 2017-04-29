/*
 * map.h
 *
 *  Created on: Dec 12, 2016
 *      Author: mufferm
 */

#ifndef MAP_H_
#define MAP_H_

class Map {
public:
	
	struct single_landmark_s{

		int id_i ; // Landmark ID
		float x; // Landmark x-position in the map (global coordinates)
		float y; // Landmark y-position in the map (global coordinates)
        
        template<typename T>
        double dist(const T& point) const {
            const double dx = point.x-x;
            const double dy = point.y-y;
            return sqrt(dx*dx+dy*dy);
        }
	};

	std::vector<single_landmark_s> landmark_list ; // List of landmarks in the map
    
    template<typename T>
    single_landmark_s nearest(const T& point) {
        single_landmark_s best = landmark_list.front();
        double best_dist = best.dist(point);
        for (const single_landmark_s &curr : landmark_list) {
            const double dist = curr.dist(point);
            if (dist < best_dist) {
                best_dist = dist;
                best      = curr;
            }
        }
        return best;
    }
};

#endif /* MAP_H_ */
