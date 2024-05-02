// ParticleContainer.cpp

#include "ParticleContainer.h"

// Constructor (if needed)
ParticleContainer::ParticleContainer() {
    // Initialize any other members if necessary
}

/**
 * @param particle: Particle that will be added to the
 */
void ParticleContainer::addParticle(const Particle& particle) {
    particles.push_back(particle);
    pairings.emplace_back(); // Initialize an empty pairing list for this particle
}

// Add a pairing between two particles
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
