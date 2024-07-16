//
// Created by kaiarenja on 16.05.24.
//

#ifndef CALCULATIONS_H
#define CALCULATIONS_H
#include "Container/ParticleContainer.h"
#include "utils/EpsilonSigma.h"
#include <map>

/**
 * @brief Calculations class
 */
class Calculations {
private:
  BaseParticleContainer &particles;
  double r_cutoff;
  double g_grav;
  bool smoothLJ;
  double r_l;
  double maxDistance = 0.0;
  double stiffness = 300;
  double avgBondLength = 2.2;

public:
  Calculations();
  Calculations(BaseParticleContainer &particles);

  /**
 * @brief sets the cutoff radius
 * @param r_cutoff is the cutoff radius

 */
  void setR_cutoff(double r_cutoff);

  /**
   * @brief sets the gravitational force
   * @param g_grav is the gravitational force
   */
  void setG_grav(double g_grav);

  /**
   * @brief setter for Deciding wether to have SmoothLJ or normal LJ
   */
  void setSmoothLJ(bool SLJ);

  void setR_L(double r_l);

  /**
   * @brief sets the stiffness
   * @param stif is the stiffness
   */
  void setStiffness(double stif);

  /**
   * @brief sets the avgBondLength
   * @param length is the avgBondLength
   */
  void setAvgBondLength(double length);

  /**
   * @brief the following function calculates the gravitational force
   *
   * therefore we use an interator and also for-loops
   * to calculate the physics behind the forces
   */
  void calculateF();

  /**
   *@brief calculate the Lennard-Jones-Force between the center particles
   *and the neighboourhood particles
   *(escpecially used for LC)
   *@param center are the particles in the center cell
   *@param other are the particles in the neighbour cell
   */
  void LCcalculateLJF(std::vector<Particle *> &center,
                      std::vector<Particle> &other,
                      const std::vector<EpsilonSigma> EAndS);
  /**
   *@brief calculate the Lennard-Jones-Force in the center particles
   *(escpecially used for LC)
   *called by Calculations::LCcalculateLJF(center, other)
   *@param center is the paramter used in the method above
   */
  void calculateLJFcenter(std::vector<Particle *> &center,
                          const std::vector<EpsilonSigma> EAndS);
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
  std::array<double, 3> calculateLJF(Particle *p1, Particle *p2, double e,
                                     double s);

  std::array<double, 3> calculateSmoothLJF(Particle *p1, Particle *p2, double e,
                                           double s);

  std::array<double, 3> decideForceMethod(Particle *p1, Particle *p2, double e,
                                          double s);

  std::vector<double> calculateHarmonicForce(Particle *p1, Particle *p2,
                                             double r0);

  double calcDistance(std::array<double, 3> x1, std::array<double, 3> x2);

  std::map<double, double>
  calculateLocalDensities(const std::vector<Particle> particles, double deltaR);

  std::vector<std::vector<double>>
  computeDistances(std::vector<Particle> particles);

  static double calculateDiffusion(std::vector<Particle> particles,
                                   std::vector<Particle> prevParticles);

  static double calculateDistanceBetweenParticles(Particle *p1, Particle *p2);
};

#endif // CALCULATIONS_H
