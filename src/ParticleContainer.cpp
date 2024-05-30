// ParticleContainer.cpp

#include "ParticleContainer.h"

#include <cmath>
#include <iostream>

#include "Calculations.h"

// Constructor (if needed)
ParticleContainer::ParticleContainer() {
  // Initialize any other members if necessary
}

/**
 * @param particle: Particle that will be added to the vector
 */
void ParticleContainer::addParticle(const Particle &particle) {
  particles.push_back(particle);
  pairings.emplace_back(); // Initialize an empty pairing list for this particle
}

/**
 * @brief add a Pairing to the 2d-vector pairings
 * @param particleIndex1: index of the first Particle in particles
 * @param particleIndex2: index of the second Particle in particles
 */
void ParticleContainer::addPairing(int particleIndex1, int particleIndex2) {
  if (particleIndex1 >= 0 && particleIndex1 < int(particles.size()) &&
      particleIndex2 >= 0 && particleIndex2 < int(particles.size())) {
    pairings.push_back(std::make_pair(particleIndex1, particleIndex2));
  }
}

std::vector<Particle> &ParticleContainer::getParticles() { return particles; }
const std::vector<Particle> &ParticleContainer::getParticles() const {
  return particles;
}

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

void ParticleContainer::resetParticles() { particles.clear(); }

void ParticleContainer::addCube(std::vector<Particle> particleCube) {
  for (size_t x = 0; x < particleCube.size(); ++x) {
    addParticle(particleCube.at(x));
  }
}

void ParticleContainer::addDisk(std::vector<Particle> particleDisk) {
  for (size_t x = 0; x < particleDisk.size(); ++x) {
    addParticle(particleDisk.at(x));
  }
}

void ParticleContainer::setParticle(Particle p, int position) {
  particles.at(position) = p;
}
std::vector<Particle> LCParticleContainer::getParticleInNeighbourhood(Cell c) {
  std::array<int, 3> id = c.getId();
  std::vector<Particle> neigbourhood;
  for(auto& c : cells) {
    int x = sqrt(std::pow(c.getId().at(0) - id.at(0),2));
    int y = sqrt(std::pow(c.getId().at(1) - id.at(1),2));
    int z = sqrt(std::pow(c.getId().at(2) - id.at(2),2 ));
    if(x <= 1 && y <= 1 && z <= 1) {
      for(auto& p: c.getParticles()) {
        neigbourhood.push_back(p);
      }
    }
  }
  return neigbourhood;
};
Cell &LCParticleContainer::getCellById(std::array<int, 3> id) {
  for(auto& c : cells) {
    if(c.getId() == id) {
      return c;
    }
  }
  Cell c({-1,-1,-1});
  return c;
}

void LCParticleContainer::realocateParticles(int handle_out_of_border) {
  //handle_out_of_border == 1 => delete Particles outside
  //handle_out_of_border == 2 => bounce back
  if(handle_out_of_border == 1) {
    //finds right cell for every particle
    std::vector<Particle> all_particles;
    for(auto& c:cells) {
      for(auto& p : c.getParticles()) {
       all_particles.push_back(p);
      }
    }
    for(auto& p : all_particles) {
      int x = floor(p.getX().at(0)/cell_size.at(0));
      int y = floor(p.getX().at(1)/cell_size.at(1));
      int z = floor(p.getX().at(2)/cell_size.at(2));
      Cell c = getCellById({x,y,z});
      std::array<int, 3> wrong_id = {-1,-1,-1};
      if(c.getId() != wrong_id) {
        //adds p to right cell in case it is not out of the boundarys
        p.setF({0.0,0.0,0.0});
        c.addParticle(p);
      }
    }
  } else {

  }
}
void LCParticleContainer::generateCells(int size_x, int size_y, int size_z, double r_cutoff) {
  int count_x = floor(size_x/r_cutoff);
  int count_y = floor(size_y/r_cutoff);
  int count_z = floor(size_z/r_cutoff);
  double cell_size_x = size_x/count_x;
  double cell_size_y = size_y/count_y;
  double cell_size_z = size_z/count_z;
  for(int x = 0; x < count_x; x++) {
    for(int y = 0; y < count_y; y++) {
      for(int z = 0; z < count_z; z++) {
        Cell c({x,y,z});
      }
    }
  }
  cell_size={cell_size_x,cell_size_y,cell_size_z};
}

void LCParticleContainer::handleLJFCalculation() {
  realocateParticles(1);
  for(auto& c:cells) {
    std::vector<Particle> neighbourhood = getParticleInNeighbourhood(c);
    Calculations calc{};
    calc.LCcalculateLJF(neighbourhood);
  }
}




