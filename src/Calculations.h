//
// Created by kaiarenja on 16.05.24.
//

#ifndef CALCULATIONS_H
#define CALCULATIONS_H
#include "Container/ParticleContainer.h"

/**
 * @brief Calculations class
 */
class Calculations {
private:
  BaseParticleContainer &particles;
  static constexpr double epsilion = 5;
  static constexpr double sigma = 1;
  double r_cutoff;
  double g_grav;

public:
  Calculations();
  Calculations(BaseParticleContainer &particles);

 void setR_cutoff(double r_cutoff);

 void setG_grav(double g_grav);

 /**
 * @brief the following function calculates the gravitational force
 *
 * therefore we use an interator and also for-loops
 * to calculate the physics behind the forces
 */
  void calculateF();
  /**
   * @brief calculate the Lennard-Jones-Force for all particles
   */
  void calculateLJF();
 /**
   *@brief calculate the Lennard-Jones-Force between the center particles
   *and the neighboourhood particles
   *(escpecially used for LC)
   *@param center are the particles in the center cell
   *@param other are the particles in the neighbour cell
   */
  void LCcalculateLJF(std::vector<Particle *> &center,
                      std::vector<Particle > &other);
 /**
   *@brief calculate the Lennard-Jones-Force in the center particles
   *(escpecially used for LC)
   *called by Calculations::LCcalculateLJF(center, other)
   *@param center is the paramter used in the method above
   */
  void calculateLJFcenter(std::vector<Particle *> &center);
 /**
* @brief the following function calculates the new positions
*
* therefore we use a for-loops
* to calculate the physics behind the positions (Velocity-Störmer-Verlet)
*/
  void calculateX(double delta_t);
 /**
* @brief the following function calculates the new velocities
*
* therefore we use a for-loops
* to calculate the physics behind the velocities (Velocity-Störmer-Verlet)
*/
  void calculateV(double delta_t);

 /**
   *@brief calculate the Lennard-Jones-Force between two particles
   *(escpecially used for LC)
   */
  std::array<double, 3> calculateLJF(Particle *p1, Particle *p2);
};

#endif // CALCULATIONS_H
