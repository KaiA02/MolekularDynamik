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
  BaseParticleContainer &particles;
  static constexpr double epsilion = 5;
  static constexpr double sigma = 1;

public:
  Calculations();
  Calculations(BaseParticleContainer &particles);

  /**
   *@brief calculate the force for all particles
   */
  void calculateF();
  /**
   * @brief calculate the Lennard-Jones-Force for all particles
   */
  void calculateLJF();
  void LCcalculateLJF(std::vector<Particle *> &center,
                      std::vector<Particle> &other);
  void calculateLJFcenter(std::vector<Particle *> &center);
  /**
   * @brief calculate the position for all particles
   */
  void calculateX(double delta_t);
  /**
   *@brief calculate the velocity for all particles
   */
  void calculateV(double delta_t);

  std::array<double, 3> calculateLJF(Particle *p1, Particle *p2);
};

#endif // CALCULATIONS_H
