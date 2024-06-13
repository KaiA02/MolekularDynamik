#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include <optional>

#include "../Particle.h"
#include <vector>

#include "../LinkedCell/Cell.h"



/**
 * @brief The ParticleContainer class
 * This class is a container for Particles. It contains a vector of Particles
 * and a 2D vector of pairings.
 */
class ParticleContainer {
public:
  ParticleContainer();
 /**
  * @param particle: Particle that will be added to the vector
  */
  void addParticle(const Particle &particle);
  void setParticle(Particle p, int position);
  void resetParticles();
  std::vector<Particle> &getParticles();
  const std::vector<Particle> &getParticles() const;
  std::vector<Cell> getCells();
  int size() const;

  std::vector<Particle>::iterator begin();
  std::vector<Particle>::iterator end();
  std::vector<Particle>::const_iterator begin() const;
  std::vector<Particle>::const_iterator end() const;

  void handleLJFCalculation();

  /**
   * @brief adds a cuboid of particles to the vector particles
   * @param particleCube cuboid generated by the ParticleGenerator
   */
  void addMultipleParticles(std::vector<Particle> particleCube);

  /**
   * @brief adds a disk of particles to the vector particles
   * @param particleDisk disk generated by the ParticleGenerator
   */
  // void addDisk(std::vector<Particle> particleDisk);

protected:
 std::vector<Particle> particles;
};

#endif // PARTICLECONTAINER_H
