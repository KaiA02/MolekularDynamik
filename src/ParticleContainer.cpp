// ParticleContainer.cpp

#include "ParticleContainer.h"

#include <iostream>

// Constructor (if needed)
ParticleContainer::ParticleContainer() {
    // Initialize any other members if necessary
}

/**
 * @param particle: Particle that will be added to the vector
 */
void ParticleContainer::addParticle(const Particle& particle) {
    particles.push_back(particle);
    pairings.emplace_back(); // Initialize an empty pairing list for this particle
}

/**
* @brief add a Pairing to the 2d-vector pairings
 * @param particleIndex1: index of the first Particle in particles
 * @param particleIndex2: index of the second Particle in particles
 */
void ParticleContainer::addPairing(int particleIndex1, int particleIndex2) {
    if (particleIndex1 >= 0 && particleIndex1 < particles.size() &&
        particleIndex2 >= 0 && particleIndex2 < particles.size()) {
        pairings[particleIndex1].push_back(particleIndex2);
        pairings[particleIndex2].push_back(particleIndex1);
        }
}

const std::vector<Particle>& ParticleContainer::getParticles() const {
        return particles;
    }

int ParticleContainer::size() const {
    return particles.size();
}

std::vector<Particle>::iterator ParticleContainer::begin() {
    return particles.begin();
}

std::vector<Particle>::iterator ParticleContainer::end() {
    return particles.end();
}

std::vector<Particle>::const_iterator ParticleContainer::begin() const {
    return particles.begin();
}

std::vector<Particle>::const_iterator ParticleContainer::end() const {
    return particles.end();
}

void ParticleContainer::addCube(std::vector<Particle> particleCube) {
    std::cout << particleCube.size() << ": particles are in the cube ------------------------------------------> " << std::endl;
    for(size_t x = 0; x < particleCube.size(); ++x) {
        addParticle(particleCube.at(x));
        std::cout << particles.size() << "particles added " << std::endl;


    }
    std::cout << particles.size() << " particles added to the container ------------------------------------------> " << std::endl;
}