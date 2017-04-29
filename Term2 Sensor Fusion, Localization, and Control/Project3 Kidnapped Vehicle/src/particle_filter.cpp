/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

// translate the observation to the particle and apply rotation matrix https://en.wikipedia.org/wiki/Rotation_matrix
LandmarkObs Particle::transform(const LandmarkObs obs) const {
    LandmarkObs res = obs;
    const double c  = cos(theta);
    const double s  = sin(theta);

    res.x = x+c*obs.x-s*obs.y;
    res.y = y+s*obs.x+c*obs.y;
    return res;
}

void ParticleFilter::init(const double x, const double y, const double theta, const double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    std::default_random_engine gen;
    std::normal_distribution<double> dx(x    , std[0]);
    std::normal_distribution<double> dy(y    , std[1]);
    std::normal_distribution<double> dt(theta, std[2]);
    for (int i=0; i<num_particles; ++i) {
        const Particle p = {i, dx(gen), dy(gen), dt(gen), 1./num_particles};
        particles.push_back(p);
    }
    
    is_initialized = true;
}

void ParticleFilter::prediction(const double delta_t, const double std_pos[], const double velocity, const double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    std::default_random_engine gen;
    std::normal_distribution<double> dx(0, std_pos[0]);
    std::normal_distribution<double> dy(0, std_pos[1]);
    std::normal_distribution<double> dt(0, std_pos[2]);
    for (Particle &p : particles) {
        if (fabs(yaw_rate) > 0.00001) {
            p.x += velocity/yaw_rate*(sin(p.theta+yaw_rate*delta_t)-sin(p.theta));
            p.y += velocity/yaw_rate*(cos(p.theta)-cos(p.theta+yaw_rate*delta_t));
        }   else {
            p.x += velocity*cos(p.theta)*delta_t;
            p.y += velocity*sin(p.theta)*delta_t;
        }
        p.x     += dx(gen);
        p.y     += dy(gen);
        p.theta += dt(gen) + yaw_rate*delta_t;
    }
}

void ParticleFilter::updateWeights(const double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

    double total_weight = 0;
    for (Particle &p : particles) {
        // the new weight of this particle is the product of best measurement weights
        p.weight = 1;
        for (const LandmarkObs &raw_obs : observations) {
            const LandmarkObs obs             = p.transform(raw_obs);
            const Map::single_landmark_s best = map_landmarks.nearest(obs);
            p.weight *= weight(obs, best, std_landmark[0], std_landmark[0]);
        }
        total_weight += p.weight;
    }

    // normalize weights
    for (Particle &p : particles) {
        p.weight /= total_weight;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::vector<double> weights;
    for (const Particle &p : particles) {
        weights.push_back(p.weight);
    }
    std::default_random_engine gen;
    std::discrete_distribution<> d(weights.begin(), weights.end());
    const auto old_particles = particles;
    for (int i=0; i<num_particles; ++i) {
        particles[i] = old_particles[d(gen)];
    }
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}

double ParticleFilter::weight(const double x1, const double y1, const double x2, const double y2, const double sigma_x, const double sigma_y) {
    const double dx = x1-x2, dy = y1-y2;
    return exp(-0.5*(dx*dx/sigma_x/sigma_x+dy*dy/sigma_y/sigma_y)) / 2/M_PI/sigma_x/sigma_y;
}
