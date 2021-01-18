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
#include <random> // Need this for sampling from distributions

#include "helper_functions.h"

using std::normal_distribution;
using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   *  
   * 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 25;  // TODO: Set the number of particles

  double std_x = std[0], std_y= std[1], std_theta = std[2]; // Standard deviations for x, y, and theta
  // intialize the noise distibution
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  // Add random noise to each particle
  std::default_random_engine gen;
  for (int i = 0; i < num_particles; ++i)
  {
    Particle particle;  // Initilize one particle
    particle.id = i;
    particle.x  = dist_x(gen);
    particle.y  = dist_y(gen);
    particle.theta  = dist_y(gen);
    particle.weight = 1.0;

    particles.push_back(particle);  // push back to the vector of particles
  }
  is_initialized = true;
}

  void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                  double velocity, double yaw_rate)
  {
    /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    std::default_random_engine gen;

    // Adding noise to the prediction
    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);

    double x_f, y_f, theta_f;

    for (int i = 0; i < num_particles; ++i)
    {
      Particle particle = particles[i];
      double x_0 = particle.x;
      double y_0 = particle.y;
      double theta_0 = particle.theta;

      if (yaw_rate < 0.1)
      {
        x_f = x_0 + (velocity) * (sin(theta_0) - sin(theta_0));
        y_f = y_0 + (velocity) * (cos(theta_0) - cos(theta_0));
        theta_f = theta_0;
      }
      else
      {
        x_f = x_0 + (velocity / yaw_rate) * (sin(theta_0 + (yaw_rate * delta_t)) - sin(theta_0));
        y_f = y_0 + (velocity / yaw_rate) * (cos(theta_0 + (yaw_rate * delta_t)) - cos(theta_0 + (yaw_rate * delta_t)));
        theta_f = yaw_rate * delta_t;
      }

      // Adding noise
      particle.x = x_f + dist_x(gen);
      particle.y = y_f + dist_y(gen);
      particle.theta = theta_f + dist_theta(gen);

      particles[i] = particle;
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
  /* 
  Ziel: predict ID = Observed machting ID.
  Steps: 1. For Loop over all predicition.
            Takin in the prediction distians to all observation. Putting this in to an Vector.
            RAnking them in a vector to the samles error.
            Associate the smales number to the pediction 
  */

 // helper for calculating the error
  double distance;
  double low_distance;


  for (int i = 0; i < observations.size(); ++i){
    low_distance = 100000000;
    for (int j = 0; j < predicted.size(); ++j){

      distance = sqrt(pow((observations[i].x - predicted[j].x), 2) + pow((observations[i].y - predicted[j].y), 2));

      if (distance < low_distance){
        low_distance = distance;
        observations[i].id = j;
      }
    }
    // Sort the vector, samllest error to the beginning
    // std::sort(distances.begin(), distances.end(), distances.weight());
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

  vector<LandmarkObs> Obs_transformed;
  LandmarkObs obs;

  for (int i=0; i<particles.size(); ++i){
    double x_part, y_part, x_obs, y_obs, theta, weight;

    x_part = particles[i].x;
    y_part = particles[i].y;
    theta = particles[i].theta;
    particles[i].weight = 0;          // Importend, so that every new call will distribute the weights new

    for (int j = 0; j < observations.size(); ++j){

      obs.x = observations[j].x;
      obs.y = observations[j].y;
     
      // transform to map x coordinate
      double x_map;
      x_map = x_part + (cos(theta) * obs.x) - (sin(theta) * obs.y);

      // transform to map y coordinate
      double y_map;
      y_map = y_part + (sin(theta) * obs.x) + (cos(theta) * obs.y);

      obs.x = x_map;
      obs.y = y_map;
      obs.id = obs.y = observations[j].id;

      Obs_transformed.push_back(obs);
    }

   // declearing Vector for landmarks in vehicle sensor range 
   vector<LandmarkObs> landmarks_range;
   LandmarkObs lm; 
   for (int k = 0; k < map_landmarks.landmark_list.size(); ++k)
   {
     double x_lm, y_lm;

     x_lm = map_landmarks.landmark_list[k].x_f;
     y_lm = map_landmarks.landmark_list[k].y_f;

     // Euklidian distanance between partikel and Landmark
     double distance = dist(x_lm, y_lm, x_part, y_part);
  
     if (distance < sensor_range){
       lm.x = x_lm;
       lm.y = y_lm;
       lm.id = map_landmarks.landmark_list[k].id_i;
       landmarks_range.push_back(lm);
     }
   }

   dataAssociation(landmarks_range, Obs_transformed);

    // Updating the Particle weight
    double sig_x = std_landmark[0], sig_y = std_landmark[1];  // Standart Deviation for x and y 
    // Calculation the weight
    double mu_x, mu_y ;
    // define outputs for observation
    weight = 1.0;
    for (int l = 0; l < landmarks_range.size(); ++l){
      mu_x = landmarks_range[l].x;
      mu_y = landmarks_range[l].y;
      for (int m = 0; m < Obs_transformed.size(); ++m)
      {
        // checking that the Observation fits to the landmark. The nears assosiation is done in dataAssociation
        if (landmarks_range[l].id == Obs_transformed[m].id){
          x_obs = Obs_transformed[m].x;
          y_obs = Obs_transformed[m].y;
        
          weight *= multiv_prob(sig_x, sig_y, x_obs, y_obs, mu_x, mu_y);
        }
      }
    }
    // Update of particle weight
    particles[i].weight = weight;
    weights.push_back(weight);
  }
  // Normalize the particle weights
   double sum_weight = std::accumulate(weights.begin(), weights.end(),0);

   for (int i = 0; i < particles.size(); ++i)
  {
    particles[i].weight /= sum_weight;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  std::discrete_distribution<int> distribution(particles.begin(), particles.end());
  std::vector<Particle> weighted(num_particles);

  for (int i = 0; i < num_particles; ++i){
    int j = distribution(gen);
    weighted.at(i) = particles.at(j);
  }

  particles = weighted;
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