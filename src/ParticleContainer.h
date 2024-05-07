
#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include <vector>
#include "Particle.h"
/**
 * @brief The ParticleContainer class
 * This class is a container for Particles. It contains a vector of Particles and a 2D vector of pairings.
 */
class ParticleContainer {
public:
    ParticleContainer();
    /**
    * @brief add a Particle to the vector
    * @param Particle particle
*/
    void addParticle(const Particle& particle);

    /*
    @brief add a Pairing to the 2d-vector pairings
    @param index1: index of the first Particle in particles index2: index of the second Particle in particles
     */
    void addPairing(int particleIndex1, int particleIndex2);
    const std::vector<Particle>& getParticles() const;
	int size() const;

	std::vector<Particle>::iterator begin();
    std::vector<Particle>::iterator end();
    std::vector<Particle>::const_iterator begin() const;
    std::vector<Particle>::const_iterator end() const;

private:
    std::vector<Particle> particles;
    std::vector<std::vector<int>> pairings; // 2D vector of pairings
};

#endif // PARTICLECONTAINER_H
