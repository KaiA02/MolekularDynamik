#include "ParticleContainer.h"
#include "../Calculations.h"
#include "../LinkedCell/Cell.h"
#include "spdlog/spdlog.h"
#include <cmath>
#include <iostream>

// Constructor for ParticleContainer
ParticleContainer::ParticleContainer() {
  // Initialize any other members if necessary
}

void ParticleContainer::addParticle(const Particle &particle) {
  particles.push_back(particle);
}

std::vector<Particle> &ParticleContainer::getParticles() { return particles; }

const std::vector<Particle> &ParticleContainer::getParticles() const {
  return particles;
}
std::vector<Cell> ParticleContainer::getCells() { return {}; }

int ParticleContainer::size() const { return particles.size(); }

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

void ParticleContainer::handleLJFCalculation() {}

void ParticleContainer::resetParticles() { particles.clear(); }

void ParticleContainer::addMultipleParticles(std::vector<Particle> particleCube) {
  for (size_t x = 0; x < particleCube.size(); ++x) {
    addParticle(particleCube.at(x));
  }
}

void ParticleContainer::setParticle(Particle p, int position) {
  particles.at(position) = p;
}








