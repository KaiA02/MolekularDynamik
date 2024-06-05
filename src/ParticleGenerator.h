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
   */
  void generateCuboid(const Particle &start, int n1, int n2, int n3,
                      double distance, double meanVelocity, int dimension);
  std::vector<Particle>& getAllParticles();

    /**
      * @brief Generate a disk of particles with the radius and a cutoff value, from a
      * center particle with a given distance
      * @param center center particle with inital velocity
      * @param radius is the radius from the center on
      * @param distance distance between the particles
      */
    void generateDisk(const Particle &center, int radius,
                        double distance, int dimension);
};

#endif // PARTICLEGENERATOR_H
