//
// Created by jh on 04.07.2024.
//

#ifndef MEMBRANE_H
#define MEMBRANE_H
#include "../Calculations.h"
#include "../Particle.h"
#include "MembranePair.h"
#include <vector>
/**
 *@brief Membrane class for generating a basic membrane
 */
class Membrane {
public:
  Membrane();

  /**
   * @brief constructor for Membrane
   * @param particles to be included in the membrane
   * @param distance the distance between the particles
   * @param fUp the force upwards which is pulling the membrane
   */
  Membrane(std::vector<Particle *> particles, double distance, double fUp);

  /**
   * @brief applies the movement of the 4 moving particles -> pulls the membrane
   * with the 4 particles
   */
  void applyMovement();

  /**
   * @brief calculates the forces between the particles
   * @param calc the calculations object
   */
  void stabilizeMembrane(Calculations &calc);

  /**
   * @brief sets the moving particles that have a force already at the beginning
   * @param movPart the particles that are moving
   */
  void setMovingParticles(std::vector<Particle *> movPart);

  /**
   *
   * @return the number of moving particles
   */
  int getMovingParticleCount();

  /**
   * @brief calculates the average size of the MembranePairs
   */
  void getAveragePairSize();

  /**
   * @brief calculates the distance between two particles
   * @param x1 position of the first particle
   * @param x2 position of the second particle
   * @return the distance between the two particles
   */
  double calcDistance(std::array<double, 3> x1, std::array<double, 3> x2);

private:
  std::vector<MembranePair> pairs;
  double distance;
  std::vector<Particle *> movingParticles;
  double forceUpwards;
};
#endif // MEMBRANE_H
