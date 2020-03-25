/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution; // Generates a Gaussian Distribution 

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 10;  // TODO: Set the number of particles
  // Create the variable for Randomness
  std::default_random_engine gen;
  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std[0]);
  // This line creates a normal (Gaussian) distribution for y
  normal_distribution<double> dist_y(y, std[1]);
  // This line creates a normal (Gaussian) distribution for theta
  normal_distribution<double> dist_theta(theta, std[2]);
  // Define a momentary particle
  Particle init_particle;
  //Initialize each particle with an uncertainty.
  for(int i; i < num_particles; i++){
    init_particle.x = dist_x(gen); 
    init_particle.y = dist_y(gen);
    init_particle.theta = dist_theta(gen);
    init_particle.weight = 1.00;
    particles.push_back(init_particle);
  }
  is_initialized = true; // Set the initialization flag
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // Create the variable for Randomness
  std::default_random_engine gen;
  // Variables to store the motion model estimated positions
  double x_f, y_f, theta_f; 
      
  for(int i=0; i < num_particles; i++){
    // According to the Bicycle motion model we have the following equations for the prediction (x,y,theta)
    x_f = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
    y_f = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
    theta_f = particles[i].theta + yaw_rate*delta_t;
    // These lines creates a normal (Gaussian) distribution for x, y, and theta
    normal_distribution<double> dist_x(x_f, std_pos[0]);
    normal_distribution<double> dist_y(y_f, std_pos[1]);
    normal_distribution<double> dist_theta(theta_f, std_pos[2]);
    //Pick a random Value from the generated Gaussians and update the position of the particle
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  //Momentary variables for the Ecluidean distance
  double x1, y1, x2, y2,eucli_distance;
  //Momentary variables to store the best score
  int best_id;
  double best_x,best_y;
  //Compare each of the predicted measurements with the observation measurements
  for(int i = 0;i < predicted.size();i++){
    x1 = predicted[i].x;
    y1 = predicted[i].y;
    double best_eucli_distance = 100;
    for(int j = 0; j < observations.size();j++){
       x2 = observations[j].x;
       y2 = observations[j].y;
       //Calculate the Euclidean distance
       eucli_distance = dist(x1, y1, x2, y2);
       // If the distance is shorter than the previous one update the new distance
       if(eucli_distance < best_eucli_distance ){
         best_eucli_distance = eucli_distance;
         best_id = observations[j].id;
         best_x = observations[j].x;
         best_y = observations[j].y;
       }
    }
    // Assign the value the observation value to the predicted value.
    predicted[i].id = best_id;
    predicted[i].x = best_x;
    predicted[i].y = best_y;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}