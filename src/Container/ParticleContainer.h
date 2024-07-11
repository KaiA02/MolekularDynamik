#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include <optional>

#include "../Particle.h"
#include <vector>

#include "BaseParticleContainer.h"
#include "../LinkedCell/Cell.h"



/**
 * @brief The ParticleContainer class
 * This class is a container for Particles. It contains a vector of Particles
 * and a 2D vector of pairings.
 */
class ParticleContainer : public BaseParticleContainer {
public:
  ParticleContainer();

 /**
  *@brief adds a Particle to the vector particles
  * @param particle: Particle that will be added to the vector
  */
  void addParticle(const Particle &particle) override;

  /**
   * @brief inserts a particle at index position of the particles vector
   * @param p the new particle
   * @param position the index where the particle will be inserted
   */
  void setParticle(Particle p, int position) override;

  /**
   * @brief resets all particles in the container by clearing the vector
   * particles
   */
  void resetParticles() override;

  /**
   *
   * @return reference to the vector particles
   */
  std::vector<Particle> &getParticles() override;
  const std::vector<Particle> &getParticles() const override;

  /**
   *
   * @return vector of cells
   */
  std::vector<Cell> getCells();

  /**
   *
   * @return the number of particles in the container
   */
  int size() const override;

  std::vector<Particle>::iterator begin() override;
  std::vector<Particle>::iterator end() override;
  std::vector<Particle>::const_iterator begin() const override;
  std::vector<Particle>::const_iterator end() const override;

  /**
   * @brief handles the calculation of the Lennard-Jones force
   */
 void handleLJFCalculation() override;

  /**
   * @brief adds a vector of particles to the container
   * @param particleCube vector of particles to be added to the container
   */
  void addMultipleParticles(std::vector<Particle> particleCube);

};



#endif // PARTICLECONTAINER_H
