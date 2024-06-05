#include "ParticleContainer.h"
#include <cmath>
#include <iostream>
#include "Calculations.h"
#include "LinkedCell/Cell.h"
#include "spdlog/spdlog.h"

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


std::vector<Particle>& ParticleContainer::getParticles() {
    return particles;
}

const std::vector<Particle> &ParticleContainer::getParticles() const {
    return particles;
}
std::vector<Cell> ParticleContainer::getCells() {
  return {};
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

void ParticleContainer::addMultipleParticles(std::vector<Particle> particleCube) {
    for (size_t x = 0; x < particleCube.size(); ++x) {
        addParticle(particleCube.at(x));
    }
}



void ParticleContainer::setParticle(Particle p, int position) {
    particles.at(position) = p;
}
std::vector<Particle> LCParticleContainer::getParticleInNeighbourhood(std::array<int, 3> id) {
  std::vector<Particle> neigbourhood;
  int particleCount=0;
  for(auto &cell : cells) {
    if(cell.getParticles().size() > 0){
      int x = cell.getId().at(0) - id.at(0);
      if(x <= 1 && x >= -1) {
        int y = cell.getId().at(1) - id.at(1);
        if(y <= 1 && y >= -1) {
          int z = cell.getId().at(2) - id.at(2);
          if(z <= 1 && z >= -1) {
            for(auto &p: cell.getParticles()) {
              neigbourhood.push_back(p);
              particleCount++;
            }
          }
        }
      }
    }
  }
  //spdlog::info("neighbourhood has {} particles", particleCount);
  if(particleCount <= 0) {
    return {};
  } else {
    spdlog::info("there was a neigbourhood of {} particles", particleCount);
    return neigbourhood;

  }

}
Cell& LCParticleContainer::getCellById(std::array<int, 3> id) {
  int failCount=0;
  if(cells.size() != 0) {
    for(auto &c : cells) {
      if(c.getId() == id) {
        spdlog::info("{} fails", failCount);
        return c;
      }
      failCount++;
    }
    spdlog::info("{} fails", failCount);
  }// Return nullptr if cell is not found
}

void LCParticleContainer::realocateParticles(int handle_out_of_border) {
  if(handle_out_of_border == 1) {
    std::vector<Particle> all_particles_in_cells;
    for(auto& c:cells) {
      for(auto& p : c.getParticles()) {
        all_particles_in_cells.push_back(p);
      }
    }
    for(auto& c:cells) {
      c.emptyCell();
    }
    for(auto& p : all_particles_in_cells) {
      int x = floor(p.getX().at(0)/cell_size.at(0));
      int y = floor(p.getX().at(1)/cell_size.at(1));
      int z = floor(p.getX().at(2)/cell_size.at(2));
      Cell c = getCellById({x,y,z});
      std::array<int, 3> id = {-1,-1,-1};
      if(c.getId() != id) {
        // adds p to the right cell in case it is not out of the boundaries
        p.setF({0.0,0.0,0.0});
        c.addParticle(p);
      }
    }
  } else {
    // Handle other cases
  }
}
void LCParticleContainer::fillCellsWithParticles() {
  for(auto p: particles) {
    addParticleToCell(p);
  }
}

void LCParticleContainer::generateCells(int size_x, int size_y, int size_z, double r_cutoff) {
  if(r_cutoff>0) {
    int count_x = floor(size_x / (r_cutoff*2));
    int count_y = floor(size_y / (r_cutoff*2));
    int count_z = floor(size_z / (r_cutoff*2));
    if(size_z < (r_cutoff*2)) {
      count_z = 1;
    }
    if(count_x > 0 && count_y > 0 && count_z > 0){
      double cell_size_x = size_x / count_x;
      double cell_size_y = size_y / count_y;
      double cell_size_z = size_z/ count_z;
      if(count_z == 1) {
        cell_size_z = r_cutoff*2;
      }
      cells.clear(); // Clear existing cells before generating new ones
      for(int x = 0; x < count_x; x++) {
        for(int y = 0; y < count_y; y++) {
          for(int z = 0; z < count_z; z++) {
            cells.push_back(Cell({x, y, z})); // Add generated cells to cells vector
          }
        }
      }
      cell_size = {cell_size_x, cell_size_y, cell_size_z};
    } else {
      spdlog::info("negative value detected in: generate Cells(count_x etc.) {}, {}, {}", count_x, count_y, count_z);
    }
  } else {
    spdlog::info("negative value detected in: generate Cells(r_cutoff) {}", r_cutoff);
  }
}

void LCParticleContainer::handleLJFCalculation() {
  realocateParticles(1);
  int empty_counter = 0;
  for(auto& c:cells) {
    if(!c.isEmpty()) {
      spdlog::info("cell is not empty------------");
      std::vector<Particle> neighbourhood = getParticleInNeighbourhood(c.getId());
      spdlog::info("calculated Neigbours");
      LCParticleContainer container;
      Calculations calc(container);
      calc.LCcalculateLJF(neighbourhood);
      spdlog::info("calculated LJF for Neigbours");
    } else {
      empty_counter ++;
    }
  }
  spdlog::info("there where {} empty cells", empty_counter);
}
void LCParticleContainer::addParticleToCell(const Particle p) {
  int x = int(floor(p.getX().at(0) / cell_size.at(0)));
  int y = int(floor(p.getX().at(1) / cell_size.at(1)));
  int z = int(floor(p.getX().at(2) / cell_size.at(2)));
  Cell cell = getCellById({x, y, z});
  std::array<int, 3> id = {-1,-1,-1};
  if (cell.getId() != id) {
    cell.addParticle(p);

  } else {
  }
}
std::vector<Particle>& LCParticleContainer::getParticles() {
  return particles;
}
void LCParticleContainer::addMultipleParticles(std::vector<Particle> newParticles) {
  for(auto p: newParticles) {
    particles.push_back(p);
    addParticleToCell(p);

  }
}
void LCParticleContainer::addParticle(Particle p) {
  particles.push_back(p);
  addParticleToCell(p);
}
std::vector<Cell> LCParticleContainer::getCells() {
  return cells;
}






