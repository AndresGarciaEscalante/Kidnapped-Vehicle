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
  num_particles = 1000;  // TODO: Set the number of particles
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
  // Define the size of the weights
  weights.resize(num_particles);
  //Initialize each particle with an uncertainty.
  for(int i; i < num_particles; i++){
    init_particle.id = i;
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
    if(yaw_rate != 0){
      x_f = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      y_f = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      theta_f = particles[i].theta + yaw_rate*delta_t;
    }
    else{
      x_f = particles[i].x + velocity*delta_t*cos(particles[i].theta);
      y_f = particles[i].y + velocity*delta_t*sin(particles[i].theta);
      theta_f = particles[i].theta;
    }
    
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
  int past_j;
  //Compare each of the predicted measurements with the observation measurements
  //std::cout << observations.size() << "   " << predicted.size() <<std::endl;
  for(int i = 0;i < observations.size();i++){
    x1 = observations[i].x;
    y1 = observations[i].y;
    double best_eucli_distance = 1000;
    for(int j = 0; j < predicted.size();j++){
      x2 = predicted[j].x;
      y2 = predicted[j].y;
      //Calculate the Euclidean distance
      eucli_distance = dist(x1, y1, x2, y2);
      // If the distance is shorter than the previous one update the new distance
      if(eucli_distance < best_eucli_distance ){
        best_eucli_distance = eucli_distance;
        past_j = j;
      }
    }
    observations[i].id = predicted[past_j].id;
    predicted.erase(predicted.begin()+past_j);
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
  // Momentary variables for the Homogeneous Transformation
  LandmarkObs temp;
  double total_weight=0.0;
  for(int i = 0; i < particles.size(); i++){
    vector<LandmarkObs> observations_map_cordi;
    //STEP 1  : Applying the homoegeneous transformation to convert the car's observations into world coordinates
    for(int j = 0; j < observations.size(); j++){
      temp.x = particles[i].x + (cos(particles[i].theta)*observations[j].x) - (sin(particles[i].theta)*observations[j].y);
      temp.y = particles[i].y + (sin(particles[i].theta)*observations[j].x) + (cos(particles[i].theta)*observations[j].y);
      observations_map_cordi.push_back(temp);   
    } 

    //STEP 2: Select the predicted landmarks that are in range of the sensor
    double x,y;
    int id; 
    vector<LandmarkObs> prediction; 
    for(int j = 0; j < map_landmarks.landmark_list.size(); j++){
      id = map_landmarks.landmark_list[j].id_i; 
      x = map_landmarks.landmark_list[j].x_f;
      y = map_landmarks.landmark_list[j].y_f;
      //std:: cout << particles[i].x << "   " << particles[i].y<< std::endl;
      if((fabs(particles[i].x - x) <= sensor_range+5) && (fabs(particles[i].y - y) <= sensor_range+5)){
        prediction.push_back(LandmarkObs {id,x,y});
      }
    }

    //STEP 3: Associate the observations_map_cordi with the prediction of landmarks 
    dataAssociation(prediction, observations_map_cordi);
    
    //STEP 4: Calculate the Particle weight
    //Momentary variables
    double multi_prob,x_cor,y_cor,mean_x,mean_y,sigma_x,sigma_y,norm,combine_prob = 1.0;
    int past_k;
    for(int j = 0; j < observations_map_cordi.size(); j++){ 
      //Search the landmark values
      for(int k = 0; k < prediction.size(); k++){
        //std::cout<< k <<std::endl;
        if(observations_map_cordi[j].id == prediction[k].id){
          past_k = k;
          break;
        }
      }
      //Assign the Multivariate-Gaussian probability density variables
      x_cor = observations_map_cordi[j].x;
      y_cor = observations_map_cordi[j].y;
      mean_x = prediction[past_k].x;
      mean_y = prediction[past_k].y;
      sigma_x = std_landmark[0];
      sigma_y = std_landmark[1];
      norm = 1/(2*M_PI*sigma_x*sigma_y);
      multi_prob = norm*exp(-( ((x_cor-mean_x)*(x_cor-mean_x))/(2*sigma_x*sigma_x) + ((y_cor-mean_y)*(y_cor-mean_y))/(2*sigma_y*sigma_y)  ));
      combine_prob*=multi_prob;
    }
    weights[i]=combine_prob; 
    total_weight+=combine_prob;
  }
  
  //STEP 5: Normalize the particle weigths 
  for(int i = 0; i < weights.size(); i++){
    weights[i] = weights[i] / total_weight;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  //Momentary variables for the Resampling wheel
  double beta = 0.0;
  vector<Particle> temp_particle;
  int index = rand() % num_particles; //Initialization of the first index position
  std::default_random_engine gen;   // Create the variable for Randomness
  double U,max_weight;
  
  for(int i = 0; i < num_particles; i++){
    //Implement the random number from 0 to 2*MaxWeigth for U
    max_weight = 2.0 * *max_element(weights.begin(), weights.end());
    U = max_weight*( (double)rand() / (double)RAND_MAX );
    beta = beta + U ;
    while(weights[index] < beta){
      beta = beta - weights[index];
      index = (index + 1) % num_particles;
    }   
    //Store all the information from the selected index 
    temp_particle.push_back(particles[index]);
  }
  //Replace the old particles with the new ones 
  for(int j = 0; j < temp_particle.size();j++){
    particles[j] = temp_particle[j];
  }
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