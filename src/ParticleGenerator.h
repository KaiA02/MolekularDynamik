//
// Created by jh on 13.05.2024.
//

#ifndef PARTICLEGENERATOR_H
#define PARTICLEGENERATOR_H
#include <vector>

#include "Particle.h"
/**
 * @brief ParticleGenerator class: used to generate a cuboid of particles
 */
class ParticleGenerator {
private:
  std::vector<Particle> allParticles;

public:
  ParticleGenerator();
  /**
   * @brief Generate a cuboid of particles with dimensions n1, n2, n3, from a
   * starting particle with a given distance and mean velocity
   * @param start starting particle
   * @param n1 x dimension of the cuboid
   * @param n2 y dimension of the cuboid
   * @param n3 z dimension of the cuboid
   * @param distance distance between the particles
   * @param meanVelocity brownian motion mean velocity
   * @param dimension defines the dimension (2D or 3D)
   * @param temp_init initial temperature
   */
  void generateCuboid(const Particle &start, int n1, int n2, int n3,
                      double distance, double meanVelocity, int dimension,
                      double temp_init);

  /**
   * @brief generates a membrane with dimensions n1, n2, n3, from a starting
   * particle. Neighbouring particles interact via the harmonic potential so
   * they have fixed neighbor relationships, but the structure of the whole
   * membrane is nevertheless flexible. making the membrane
   * @param start the start particle
   * @param n1 x dimension of the membrane
   * @param n2 y dimension of the membrane
   * @param n3 z dimension of the membrane
   * @param distance distance between the particles
   * @param meanVelocity brownian motion mean velocity
   * @param dimension defines the dimension (2D or 3D)
   * @param temp_init initial temperature
   */
  void generateMembrane(const Particle &start, int n1, int n2, int n3,
                        double distance, double meanVelocity, int dimension,
                        double temp_init);

  std::vector<Particle> &getAllParticles();

  /**
   * @brief Generate a disk of particles with the radius and a cutoff value,
   * from a center particle with a given distance
   * @param center center particle with inital velocity
   * @param radius is the radius from the center on
   * @param distance distance between the particles
   * @param meanVelocity brownian motion mean velocity
   * @param dimension defines the dimension (2D or 3D)
   * @param temp_init initial temperature
   */
  void generateDisk(const Particle &center, int radius, double distance,
                    double meanVelocity, int dimension, double temp_init);
};

#endif // PARTICLEGENERATOR_H
