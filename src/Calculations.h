//
// Created by kaiarenja on 16.05.24.
//

#ifndef CALCULATIONS_H
#define CALCULATIONS_H
#include "ParticleContainer.h"
/**
 * @brief Calculations class
 */
class Calculations {
private:
  ParticleContainer& particles;
  static constexpr double epsilion = 5;
  static constexpr double sigma = 1;

public:
  Calculations();
  Calculations(ParticleContainer &other);
  /**
   *@brief calculate the force for all particles
   */
  void calculateF();
  /**
   * @brief calculate the Lennard-Jones-Force for all particles
   */
  void calculateLJF();
  void LCcalculateLJF(std::vector<Particle>& other);
  /**
   * @brief calculate the position for all particles
   */
  void calculateX(double delta_t);
  /**
   *@brief calculate the velocity for all particles
   */
  void calculateV(double delta_t);

  ParticleContainer &getParticles();
};

#endif // CALCULATIONS_H
