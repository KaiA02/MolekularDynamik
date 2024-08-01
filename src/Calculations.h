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
  int ParallelStrategy = 0;

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

  void setParallelStrategy(int strat);

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

  /**
   * calculates the smoothed Lennard-Jones-Force between two particles
   * @param p1 the first particle
   * @param p2 the second particle
   * @param e the epsilon value
   * @param s the sigma value
   * @return the smoothed Lennard-Jones-Force between two particles
   */
  std::array<double, 3> calculateSmoothLJF(Particle *p1, Particle *p2, double e,
                                           double s);

  /**
   * @brief decides which force method to use and returns the force between two
   * particles
   * @param p1 the first particle
   * @param p2 the second particle
   * @param e the epsilon value
   * @param s the sigma value
   * @return the force between two particles
   */
  std::array<double, 3> decideForceMethod(Particle *p1, Particle *p2, double e,
                                          double s);

  /**
   * @brief calculates the harmonic force between two particles which is used
   * for the membrane
   * @param p1 the first particle
   * @param p2 the second particle
   * @param r0 average bond length
   * @return the harmonic force between two particles
   */
  std::vector<double> calculateHarmonicForce(Particle *p1, Particle *p2,
                                             double r0);

  /**
   *@brief calculates the distance between two particles
   * @param x1 the position of the first particle
   * @param x2 the position of the second particle
   * @return the distance between the two particles
   */
  double calcDistance(std::array<double, 3> x1, std::array<double, 3> x2);

  /**
   * @brief calculates the local densities of the particles with respect to the
   * interval size deltaR, for the RDF plot
   * @param particles the particles for which the local densities are calculated
   * @param deltaR the interval size
   * @return a map with the interval mapped to the respective local density
   */
  std::map<double, double>
  calculateLocalDensities(const std::vector<Particle> particles, double deltaR);

  /**
   * @brief calculates the distances between all particle pairs
   * @param particles the vector of particles
   * @return a 2D vector with the distances between all particle pairs
   */
  std::vector<std::vector<double>>
  computeDistances(std::vector<Particle> particles);

  /**
   * @brief calculates the diffusion of the particles with respect to the
   * previous particles
   * @param particles the vector of particles
   * @param prevParticles the vector of previous particles
   * @return the diffusion of the particles
   */
  static double calculateDiffusion(std::vector<Particle> particles,
                                   std::vector<Particle> prevParticles);

  /**
   * @brief calculates the distance between two particles
   * @param p1 the first particle
   * @param p2 the second particle
   * @return the distance between the two particles
   */
  static double calculateDistanceBetweenParticles(Particle *p1, Particle *p2);
};

#endif // CALCULATIONS_H
