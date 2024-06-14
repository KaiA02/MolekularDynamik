#include "ParticleContainer.h"
#include "Calculations.h"
#include "LinkedCell/Cell.h"
#include "spdlog/spdlog.h"
#include <cmath>
#include <iostream>
#include <omp.h>

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

void LCParticleContainer::setR_cutoff(double radius) {
  r_cutoff = radius;
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

Cell* LCParticleContainer::getCellById(std::array<int, 3> id) {
  if (cells.size() != 0) {
    for (auto& c : cells) {
      if (c.getId() == id) {
        return &c;
      }
    }
  }
  return nullptr;
}

void LCParticleContainer::realocateParticles() {
  //spdlog::info("we have {} cells and {} particles", cells.size(),getParticles().size());
  for (auto &c : cells) {
    c.emptyCell();
  }
  //std::vector<Particle> updatedParticles;
  for (auto &p : particles) {
    int x = floor(p.getX().at(0) / cell_size.at(0));
    int y = floor(p.getX().at(1) / cell_size.at(1));
    int z = floor(std::abs(p.getX().at(2)) / cell_size.at(2));
    p.setF({0.0, 0.0, 0.0});
    if (cellExists({x, y, z})) {
      getCellById({x,y,z})->addParticle(&p);
      //updatedParticles.push_back(p);
    } else { //Parking of deleted Particles
        p.setF({0.0,0.0,0.0});
        p.setX({-10.0,-10.0,-10.0}); //todo: adjust this Parking Spot regarding to DomainSize
        p.setV({0.0,0.0,0.0});
        p.setType(p.getType()+10);
      }
    }
    //particles = updatedParticles;
  }

void LCParticleContainer::fillCellsWithParticles() {
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
      for (int x = -1; x < cell_count[0] +1; x++) {
        for (int y = -1; y < cell_count[1] +1; y++) {
          for (int z = -1; z < cell_count[2] +1; z++) {
            bool isHalo;
            if(x == -1 || y== -1 || z == -1 || x == cell_count[0] || y == cell_count[1] || z == cell_count[2]) {
              isHalo = true;
              Cell c({x,y,z}, isHalo);
              std::array<int,3> oponentID = findOponentCellID(c.getId());
              Cell oponent(oponentID, isHalo), *o;
              c.setOposition(o);
              cells.push_back(c);
            } else {
              isHalo = false;
            }
            cells.push_back(Cell({x, y, z}, isHalo)); // Add generated cells to cells vector
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
  realocateParticles();
  int empty_counter = 0;
  for (auto &c : cells) {
    if (!c.isEmpty()) {
      std::vector<Particle> neighbourhood =
          getParticleInNeighbourhood(c.getId());
      LCParticleContainer container;
      Calculations calc(container);
      calc.setR_cutoff(r_cutoff);
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
    getCellById({x, y, z})->addParticle(&p);
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
  for(int x = -1; x < cell_count[0] +1; x++) {
    for(int y = -1; y < cell_count[1]+1; y++) {
      for(int z = -1; z < cell_count[2]+1; z++) {
        if(x<=0 || x>=(cell_count[0]-1) || y<=0 || y>=(cell_count[1]-1) || z<=0 || z>=(cell_count[2]-1)) {
          if (cellExists({x,y,z})) {
            for(auto p : getCellById({x,y,z})->getParticles()) {
              result.push_back(p);
            }
          }
        }
      }
    }
  }
  return result;
}

void LCParticleContainer::handleBoundaryAction() {
  std::vector<Particle*> boundaryparticles = getBoundaryParticles();
    for(auto p: boundaryparticles) {
      std::array<double, 6> bounds = getInfluencingBoundarysWithDistance(p);
        for(int i = 0; i < 6; i++) {
          if(boundary_types[i] == 1){} //outflow
          else if(boundary_types[i] == 2){
            // generate halo Particle
            if(bounds.at(i) != -1.0 && abs(bounds.at(i)) <= r_cutoff){ //&& abs(bounds.at(i)) <= r_cutoff
              if( i == 0 ){ //Boundary to YZ Plane
                std::array<double, 3> x_arg = {-bounds.at(i), p->getX().at(1), p->getX().at(2)};
                std::array<double, 3> v_arg = {-1*(p->getV().at(0)), p->getV().at(1), p->getV().at(2)};
                calcWithHalo(p, x_arg, v_arg);
              } else if ( i == 1 ) { //Boundary to other YZ Plane
                std::array<double, 3> x_arg = {cell_size.at(0)*cell_count.at(0) + bounds.at(i), p->getX().at(1), p->getX().at(2)};
                std::array<double, 3> v_arg = {-1*(p->getV().at(0)), p->getV().at(1), p->getV().at(2)};
                calcWithHalo(p, x_arg, v_arg);
              } else if ( i == 2 ) { //Boundary to XZ Plane
                std::array<double, 3> x_arg = {p->getX().at(0), -bounds.at(1), p->getX().at(2)};
                std::array<double, 3> v_arg = {p->getV().at(0), -1*(p->getV().at(1)), p->getV().at(2)};
                calcWithHalo(p, x_arg, v_arg);
              } else if ( i == 3 ) { //Boundary to other XZ Plane
                std::array<double, 3> x_arg = {p->getX().at(0), cell_size.at(1)*cell_count.at(1) + bounds.at(i), p->getX().at(2)};
                std::array<double, 3> v_arg = {p->getV().at(0), -1*(p->getV().at(1)), p->getV().at(2)};
                calcWithHalo(p, x_arg, v_arg);
              } else if ( i == 4 ) { //Boundary to XY Plane
                std::array<double, 3> x_arg = {p->getX().at(0), p->getX().at(1), -bounds.at(2)};
                std::array<double, 3> v_arg = {p->getV().at(0), p->getV().at(1), -1*(p->getV().at(2))};
                calcWithHalo(p, x_arg, v_arg);
              } else { // i == 5      //Boundary to other XY Plane
                std::array<double, 3> x_arg = {p->getX().at(0),p->getX().at(1), cell_size.at(2)*cell_count.at(2) + bounds.at(i)};
                std::array<double, 3> v_arg = {p->getV().at(0), -p->getV().at(1), -1*(p->getV().at(2))};
                calcWithHalo(p, x_arg, v_arg);
              }
              //Particle has same absolute velocity but
              //eg. boundary is XY, then velocity in x and y are same but in z is -z;
              //Particle has same x and y position but z is same distance from boundary but reversed
              //Particle has no force

            } // if(bounds.at(i) != -1.0 && bounds.at(i) <= (pow(2, 1/6)/2)){

          } //reflectiv
          else if(boundary_types[i] == 3) {

            //creates new Particle in OponentCell with same Atributes exept for X_Arg
            std::array<double,3> x = findOponentXYZ(p->getX());
            Particle oponent(x, p->getV(), p->getM(), p->getType());
            particles.push_back(oponent);
          } //periodic

          }// if(boundary_types[i] == 2){

      }// for(auto p: boundaryparticles) {

  }


std::array<double, 6> LCParticleContainer::getInfluencingBoundarysWithDistance(Particle * p) {
  std::array<double, 6> result = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  //get position and decide which boundary is affecting the Particle in the future
  double x = p->getX().at(0);
  double y = p->getX().at(1);
  double z = p->getX().at(2);
  if(x <= cell_size.at(0)){ //Particle close to boundary 1
    double distance1 = x;
    result.at(0) = distance1;
  }
  if(x >= (cell_size.at(0)*(cell_count.at(0)-1))) { //Particle close to boundary 2
    double distance2 = cell_size.at(0) - (x - (cell_size.at(0)*(cell_count.at(0)-1)));
    result.at(1) = distance2;
  }
  if(y <= cell_size.at(1)){ //Particle close to boundary 3
    double distance3 = y;
    result.at(2) = distance3;
  }
  if(y >= (cell_size.at(1)*(cell_count.at(1)-1))) { //Particle close to boundary 4
    double distance4 = cell_size.at(1) - (y - (cell_size.at(1)*(cell_count.at(1)-1)));
    result.at(3) = distance4;
  }
  if(z <= cell_size.at(2)){ //Particle close to boundary 5
    double distance5 = z;
    result.at(4) = distance5;
  }
  if(z >= (cell_size.at(2)*(cell_count.at(2)-1))) { //Particle close to boundary 6
    double distance6 = cell_size.at(2) - (z - (cell_size.at(2)*(cell_count.at(2)-1)));
    result.at(5) = distance6;
  }
  return result;
}

std::array<int, 3> LCParticleContainer::findOponentCellID(std::array<int, 3> id) {
  int x;
  int y;
  int z;
  if(id.at(0) == -1) {
    x = cell_count.at(0);
  }else if(id.at(0) == cell_count.at(0)) {
    x = -1;
  } else {
    x = cell_count.at(0)-id.at(0);
  }
  if(id.at(1) == -1) {
    y = cell_count.at(1);
  }else if(id.at(1) == cell_count.at(1)) {
    x = -1;
  } else {
    y = cell_count.at(1)-id.at(1);
  }
  if(id.at(2) == -1) {
    z = cell_count.at(2);
  }else if(id.at(2) == cell_count.at(2)) {
    z = -1;
  } else {
    z = cell_count.at(2)-id.at(2);
  }
  return {x, y, z};

}

std::array<double, 3> LCParticleContainer::findOponentXYZ(std::array<double, 3> XYZ) {
  double x;
  double y;
  double z;
  if(XYZ.at(0) == -1) {
    x = cell_count.at(0) * cell_size.at(0) + XYZ.at(0);
  }else if(XYZ.at(0) == cell_count.at(0)) {
    x = -1 + XYZ.at(0);
  } else {
    x = cell_count.at(0) * cell_size.at(0) - XYZ.at(0);
  }
  if(XYZ.at(1) == -1) {
    y = cell_count.at(1) * cell_size.at(1) + XYZ.at(1);
  }else if(XYZ.at(1) == cell_count.at(1)) {
    x = -1 + XYZ.at(1);
  } else {
    y = cell_count.at(1)*cell_size.at(1)-XYZ.at(1) ;
  }
  if(XYZ.at(2) == -1) {
    z = cell_count.at(2) * cell_size.at(2) + XYZ.at(2);
  }else if(XYZ.at(2) == cell_count.at(2)) {
    z = -1 + XYZ.at(2);
  } else {
    z = cell_count.at(2)*cell_size.at(2)-XYZ.at(2);
  }
  return {x, y, z};
}

void LCParticleContainer::calcWithHalo(Particle *p, std::array<double, 3> x_arg, std::array<double, 3> v_arg) {
  double m_arg = p->getM();
  int type_arg = p->getType();
  Particle haloParticle = Particle(x_arg, v_arg, m_arg, type_arg);
  Calculations calc(*this);
  calc.calculateLJF(p, &haloParticle);
}






