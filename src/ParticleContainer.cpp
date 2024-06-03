#include "ParticleContainer.h"
#include <cmath>
#include <iostream>
#include "Calculations.h"
#include "LinkedCell/Cell.h"

// Constructor for ParticleContainer
ParticleContainer::ParticleContainer() {
    // Initialize any other members if necessary
}

/**
 * @param particle: Particle that will be added to the vector
 */
void ParticleContainer::addParticle(const Particle &particle) {
    particles.push_back(particle);
}


std::vector<Particle> &ParticleContainer::getParticles() {
    return particles;
}

const std::vector<Particle> &ParticleContainer::getParticles() const {
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

void ParticleContainer::handleLJFCalculation() {
}

void ParticleContainer::resetParticles() {
    particles.clear();
}

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
std::vector<Particle> LCParticleContainer::getParticleInNeighbourhood(Cell cell) {
  std::array<int, 3> id = cell.getId();
  std::vector<Particle> neigbourhood;
  for(auto& cell : cells) {
    int x = sqrt(std::pow(cell.getId().at(0) - id.at(0), 2));
    int y = sqrt(std::pow(cell.getId().at(1) - id.at(1), 2));
    int z = sqrt(std::pow(cell.getId().at(2) - id.at(2), 2));
    if(x <= 1 && y <= 1 && z <= 1) {
      for(auto& p: cell.getParticles()) {
        neigbourhood.push_back(p);
      }
    }
  }
  return neigbourhood;
}
Cell* LCParticleContainer::getCellById(std::array<int, 3> id) {
  if(cells.size() != 0) {
    for(auto& c : cells) {
      if(c.getId() == id) {
        return &c;
      }
    }
  }
  return nullptr; // Return nullptr if cell is not found
}

void LCParticleContainer::realocateParticles(int handle_out_of_border) {
  if(handle_out_of_border == 1) {
    std::vector<Particle> all_particles;
    for(auto& c:cells) {
      for(auto& p : c.getParticles()) {
        all_particles.push_back(p);
      }
    }
    for(auto& c:cells) {
      c.emptyCell();
    }
    for(auto& p : all_particles) {
      int x = floor(p.getX().at(0)/cell_size.at(0));
      int y = floor(p.getX().at(1)/cell_size.at(1));
      int z = floor(p.getX().at(2)/cell_size.at(2));
      Cell* c = getCellById({x,y,z});
      if(c != nullptr) {
        // adds p to the right cell in case it is not out of the boundaries
        p.setF({0.0,0.0,0.0});
        c->addParticle(p);
      }
    }
  } else {
    // Handle other cases
  }
}
void LCParticleContainer::generateCells(int size_x, int size_y, int size_z, double r_cutoff) {
  int count_x = floor(size_x / r_cutoff);
  int count_y = floor(size_y / r_cutoff);
  int count_z = floor(size_z / r_cutoff);
  double cell_size_x = size_x / count_x;
  double cell_size_y = size_y / count_y;
  double cell_size_z = size_z / count_z;
  cells.clear(); // Clear existing cells before generating new ones
  for(int x = 0; x < count_x; x++) {
    for(int y = 0; y < count_y; y++) {
      for(int z = 0; z < count_z; z++) {
        cells.push_back(Cell({x, y, z})); // Add generated cells to cells vector
      }
    }
  }
  cell_size = {cell_size_x, cell_size_y, cell_size_z};
}

void LCParticleContainer::handleLJFCalculation() {
  realocateParticles(1);
  for(auto& c:cells) {
    std::vector<Particle> neighbourhood = getParticleInNeighbourhood(c);
    LCParticleContainer container;
    Calculations calc(container);
    calc.LCcalculateLJF(neighbourhood);
  }
}
void LCParticleContainer::addParticle(Particle p) {
  int x = int(floor(p.getX().at(0) / cell_size.at(0)));
  int y = int(floor(p.getX().at(1) / cell_size.at(1)));
  int z = int(floor(p.getX().at(2) / cell_size.at(2)));
  Cell* cell = getCellById({x, y, z});
  if (cell != nullptr) {
    cell->addParticle(p);
  } else {
    // Handle the case where the cell is out of bounds
    //std::cerr << "Particle out of bounds: " << p << std::endl;
  }
}





