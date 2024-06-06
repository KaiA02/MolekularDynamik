#include "ParticleContainer.h"
#include "Calculations.h"
#include "LinkedCell/Cell.h"
#include "spdlog/spdlog.h"
#include <cmath>
#include <iostream>

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

void LCParticleContainer::setBoundarys(std::array<int, 6> in) {
  boundary_types = in;
}

void ParticleContainer::setParticle(Particle p, int position) {
  particles.at(position) = p;
}

std::vector<Particle>
LCParticleContainer::getParticleInNeighbourhood(std::array<int, 3> id) {
  std::vector<Particle> neigbourhood;
  int particleCount = 0;
  for (auto &cell : cells) {
    if (cell.getParticles().size() > 0) {
      int x = cell.getId().at(0) - id.at(0);
      if (x <= 1 && x >= -1) {
        int y = cell.getId().at(1) - id.at(1);
        if (y <= 1 && y >= -1) {
          int z = cell.getId().at(2) - id.at(2);
          if (z <= 1 && z >= -1) {
            if (x != 0 || y != 0 || z != 0) {
              for (auto p : cell.getParticles()) {
                neigbourhood.push_back(*p);
                particleCount++;
              }
            }
          }
        }
      }
    }
  }
  // spdlog::info("neighbourhood has {} particles", particleCount);
  if (particleCount <= 0) {
    return {};
  } else {
    //spdlog::info("there was a neigbourhood of {} particles", particleCount);
    return neigbourhood;
  }
}

Cell &LCParticleContainer::getCellById(std::array<int, 3> id) {
  if (cells.size() != 0) {
    for (auto &c : cells) {
      if (c.getId() == id) {
        return c;
      }
    }
  }
}

void LCParticleContainer::realocateParticles() {
  //spdlog::info("we have {} cells and {} particles", cells.size(),getParticles().size());
  for (auto &c : cells) {
    c.emptyCell();
  }
  std::vector<Particle> updatedParticles;
  for (auto &p : particles) {
    int x = floor(p.getX().at(0) / cell_size.at(0));
    int y = floor(p.getX().at(1) / cell_size.at(1));
    int z = floor(std::abs(p.getX().at(2)) / cell_size.at(2));
    p.setF({0.0, 0.0, 0.0});
    if (cellExists({x, y, z})) {
      getCellById({x, y, z}).addParticle(&p);
      updatedParticles.push_back(p);
    } else {
        //when exiting at boundary 1, do what boundary 1 wants
          //case outflow: do not add Particle to List
          //case reflecting: do reflecting stuff
      }
    }
    particles = updatedParticles;
  }

void LCParticleContainer::fillCellsWithParticles() {
  int counter = 0;
  for (auto p : particles) {
    addParticleToCell(p);
  }
}

void LCParticleContainer::generateCells(int size_x, int size_y, int size_z, double r_cutoff) {
  if (r_cutoff > 0) {
    cell_count[0] = floor(size_x / (r_cutoff));
    cell_count[1] = floor(size_y / (r_cutoff));
    cell_count[2] = 1;
    if (size_z > r_cutoff) {
      cell_count[2] = floor(size_z / (r_cutoff));
    }
    if (cell_count[0] > 0 && cell_count[1] > 0 && cell_count[2] > 0) {
      double cell_size_x = size_x / cell_count[0];
      double cell_size_y = size_y / cell_count[1];
      double cell_size_z = size_z / cell_count[2];
      if (cell_count[2] == 1) {
        cell_size_z = r_cutoff;
      }
      spdlog::info("Cellcounts: {} {} {}", cell_count[0], cell_count[1], cell_count[2]);
      spdlog::info("Cellsize: {} {} {}", cell_size_x, cell_size_y, cell_size_z);
      spdlog::info("Domainsize: {} {} {}", size_x, size_y, size_z);
      spdlog::info("r_cutoff: {}", r_cutoff);
      cells.clear(); // Clear existing cells before generating new ones
      for (int x = 0; x < cell_count[0]; x++) {
        for (int y = 0; y < cell_count[1]; y++) {
          for (int z = 0; z < cell_count[2]; z++) {
            cells.push_back(
                Cell({x, y, z})); // Add generated cells to cells vector
          }
        }
      }
      cell_size = {cell_size_x, cell_size_y, cell_size_z};
    } else {
      spdlog::info(
          "negative value detected in: generate Cells(count_x etc.) {}, {}, {}",
          cell_count[0], cell_count[1], cell_count[2]);
    }
  } else {
    spdlog::info("negative value detected in: generate Cells(r_cutoff) {}",
                 r_cutoff);
  }
}

void LCParticleContainer::handleLJFCalculation() {
  std::array<int, 6> boundarys = {1,1,1,1,1,1};
  realocateParticles();
  int empty_counter = 0;
  for (auto &c : cells) {
    if (!c.isEmpty()) {
      std::vector<Particle> neighbourhood =
          getParticleInNeighbourhood(c.getId());
      LCParticleContainer container;
      Calculations calc(container);
      calc.LCcalculateLJF(c.getParticles(), neighbourhood);

    } else {
      empty_counter++;
    }
  }
  spdlog::debug("there where {} empty cells", empty_counter);
}

bool LCParticleContainer::addParticleToCell(Particle &p) {
  int x = int(floor(p.getX().at(0) / cell_size.at(0)));
  int y = int(floor(p.getX().at(1) / cell_size.at(1)));
  int z = int(floor(p.getX().at(2) / cell_size.at(2)));
  if (cellExists({x, y, z})) {
    getCellById({x, y, z}).addParticle(&p);
    return true;
  } else {
    return false;
  }
}

std::vector<Particle> &LCParticleContainer::getParticles() { return particles; }

void LCParticleContainer::addMultipleParticles(std::vector<Particle> &newParticles) {
  for (auto &p : newParticles) {
    particles.push_back(p);
  }
}

void LCParticleContainer::addParticle(Particle p) {
  particles.push_back(p);
  addParticleToCell(p);
}

std::vector<Cell> LCParticleContainer::getCells() { return cells; }

bool LCParticleContainer::cellExists(std::array<int, 3> id) {
  for (auto c : cells) {
    if (c.getId().at(0) == id.at(0) && c.getId().at(1) == id.at(1) &&
        c.getId().at(2) == id.at(2)) {
      return true;
    }
  }
  return false;
};

void LCParticleContainer::countParticlesInCells() {
  int counter = 0;
  for(auto c: cells) {
    for(auto p: c.getParticles()) {
      counter ++;
    }
  }
  spdlog::info("there are {} particles in all the cells", counter);
} //just for debugging
std::vector<Particle *> LCParticleContainer::getBoundaryParticles() {
  std::vector<Particle*> result;
  for(int x = 0; x < cell_count[0]; x++) {
    for(int y = 0; y < cell_count[1]; y++) {
      for(int z = 0; z < cell_count[2]; z++) {
        if(x==0 || x==cell_count[0]-1 || y==0 || y==cell_count[1]-1 || z==0 || z==cell_count[2]-1) {
          for(auto p : getCellById({x,y,z}).getParticles()) {
            result.push_back(p);
          }
        }
      }
    }
  }
  return result;
}

void LCParticleContainer::handleBoundaryAction() {
  std::vector<Particle*> boundaryparticles = getBoundaryParticles();
  if(boundaryparticles.size()>0) {
    for(auto p: boundaryparticles) {
      std::vector<int> bounds = getInfluencingBoundarys(p);
      if(bounds.size() > 0) {
        // generate halo Particle
          //Particle has same absolute velocity but
            //eg. boundary is XY, then velocity in x and y are same but in z is -z;
          //Particle has same x and y position but z is same distance from boundary but reversed
          //Particle has no force

      }
    }
  }
}

std::vector<int> LCParticleContainer::getInfluencingBoundarys(Particle * p) {
  std::vector<int> result;
  //get position and decide which boundary is affecting the Particle in the future
  double x = p->getX().at(0);
  double y = p->getX().at(1);
  double z = p->getX().at(2);
  if(x <= cell_size.at(0)){ //Particle close to boundary 1
    result.push_back(1);
  }
  if(x >= (cell_size.at(0)+(cell_count.at(0)-1))) { //Particle close to boundary 2
    result.push_back(2);
  }
  if(y <= cell_size.at(1)){ //Particle close to boundary 3
    result.push_back(3);
  }
  if(y >= (cell_size.at(1)+(cell_count.at(1)-1))) { //Particle close to boundary 4
    result.push_back(4);
  }
  if(z <= cell_size.at(2)){ //Particle close to boundary 5
    result.push_back(5);
  }
  if(z >= (cell_size.at(2)+(cell_count.at(2)-1))) { //Particle close to boundary 6
    result.push_back(6);
  }
  return result;
}

