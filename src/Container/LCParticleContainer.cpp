#include "LCParticleContainer.h"
#include <cmath>
#include "spdlog/spdlog.h"
#include "../Calculations.h"

LCParticleContainer::LCParticleContainer() {
}

void LCParticleContainer::setBoundarys(std::array<int, 6> in) {
    boundary_types = in;
}
double LCParticleContainer::getR_cutoff() {
  return r_cutoff;
}


std::vector<Particle*>
LCParticleContainer::getParticleInNeighbourhood(std::array<int, 3> id) {
  std::vector<Particle*> neigbourhood;
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
                neigbourhood.push_back(p);
                particleCount++;
              }
            }
          }
        }
      }
    }
  }
  if (particleCount <= 0) {
    return {};
  } else {
    return neigbourhood;
  }
}

Cell* LCParticleContainer::getCellById(std::array<int, 3> id) {
  if (cells.size() != 0) {
    for (auto &c : cells) {
      if (c.getId() == id) {
        return &c;
      }
    }
  }
  return nullptr;
}


void LCParticleContainer::realocateParticles() {
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
      updatedParticles.push_back(p);
      getCellById({x, y, z})->addParticle(&updatedParticles.back());
    } else {
        //when exiting at boundary 1, do what boundary 1 wants
          //case outflow: do not add Particle to List
          //case reflecting: do reflecting stuff
      }
    }
    particles = updatedParticles;
  }

void LCParticleContainer::fillCellsWithParticles() {
  for (auto p : particles) {
    addParticleToCell(p);
  }
}

void LCParticleContainer::generateCells(int size_x, int size_y, int size_z, double r_cutoff_input) {
  r_cutoff = r_cutoff_input;
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
    }
  }
}

void LCParticleContainer::handleLJFCalculation() {
  realocateParticles();
  int empty_counter = 0;
  for (auto &c : cells) {
    if (!c.isEmpty()) {
      std::vector<Particle*> neighbourhood =
          getParticleInNeighbourhood(c.getId());
      Calculations calc(this);
      calc.LCcalculateLJF(c.getParticles(), neighbourhood, r_cutoff);

    } else {
      empty_counter++;
    }
  }
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
} //just for debugging
std::vector<Particle *> LCParticleContainer::getBoundaryParticles() {
  std::vector<Particle*> result;
  for(int x = 0; x < cell_count[0]; x++) {
    for(int y = 0; y < cell_count[1]; y++) {
      //for(int z = 0; z < cell_count[2]; z++) {
        if(x==0 || x==(cell_count[0]-1) || y==0 || y==(cell_count[1]-1) ) { // || z==0 || z==(cell_count[2]-1)
          for(auto p : getCellById({x,y,0})->getParticles()) {
            result.push_back(p);
          }
        }
      }
  }
  return result;
}

void LCParticleContainer::handleBoundaryAction() {
  Calculations calc(this);
  std::vector<Particle*> boundaryparticles = getBoundaryParticles();
    for(auto p: boundaryparticles) {
      std::array<double, 6> bounds = getInfluencingBoundarysWithDistance(p);
        for(int i = 0; i < 6; i++) {
          if(boundary_types[i] == 2){
            // generate halo Particle
            if(bounds.at(i) != -1.0 && bounds.at(i) <= (pow(2, 1/6)/2)){
              if( i == 0 ){ //Boundary to YZ Plane
                std::array<double, 3> x_arg = {-bounds.at(i), p->getX().at(1), p->getX().at(2)};
                std::array<double, 3> v_arg = {-1*(p->getV().at(0)), p->getV().at(1), p->getV().at(2)};
                double m_arg = p->getM();
                int type_arg = p->getType();
                Particle haloParticle = Particle(x_arg, v_arg, m_arg, type_arg);
                std::array<double, 3> force = calc.calculateLJF(&haloParticle, p, r_cutoff);
                std::array<double, 3> addedForce;
                for(int k = 0; k < 3; k++) {
                  addedForce[k] = p->getF()[k] - force[k];
                }
                p->setF(addedForce);

              } else if ( i == 1 ) { //Boundary to other YZ Plane
                std::array<double, 3> x_arg = {cell_size.at(0)*cell_count.at(0) + bounds.at(i), p->getX().at(1), p->getX().at(2)};
                std::array<double, 3> v_arg = {-1*(p->getV().at(0)), p->getV().at(1), p->getV().at(2)};
                double m_arg = p->getM();
                int type_arg = p->getType();

                Particle haloParticle = Particle(x_arg, v_arg, m_arg, type_arg);
                std::array<double, 3> force = calc.calculateLJF(&haloParticle, p , r_cutoff);
                std::array<double, 3> addedForce;
                for(int k = 0; k < 3; k++) {
                  addedForce[k] = p->getF()[k] + force[k];
                }
                p->setF(addedForce);
              } else if ( i == 2 ) { //Boundary to XZ Plane
                std::array<double, 3> x_arg = {p->getX().at(0), -bounds.at(1), p->getX().at(2)};
                std::array<double, 3> v_arg = {p->getV().at(0), -1*(p->getV().at(1)), p->getV().at(2)};
                double m_arg = p->getM();
                int type_arg = p->getType();
                Particle haloParticle = Particle(x_arg, v_arg, m_arg, type_arg);
                std::array<double, 3> force = calc.calculateLJF(&haloParticle, p, r_cutoff);
                std::array<double, 3> addedForce;
                for(int k = 0; k < 3; k++) {
                  addedForce[k] = p->getF()[k] + force[k];
                }
                p->setF(addedForce);
              } else if ( i == 3 ) { //Boundary to other XZ Plane
                std::array<double, 3> x_arg = {p->getX().at(0), cell_size.at(1)*cell_count.at(1) + bounds.at(i), p->getX().at(2)};
                std::array<double, 3> v_arg = {p->getV().at(0), -1*(p->getV().at(1)), p->getV().at(2)};
                double m_arg = p->getM();
                int type_arg = p->getType();

                Particle haloParticle = Particle(x_arg, v_arg, m_arg, type_arg);
                std::array<double, 3> force = calc.calculateLJF(&haloParticle, p, r_cutoff);
                std::array<double, 3> addedForce;
                for(int k = 0; k < 3; k++) {
                  addedForce[k] = p->getF()[k] + force[k];
                }
                p->setF(addedForce);
              } else if ( i == 4 ) { //Boundary to XY Plane
                std::array<double, 3> x_arg = {p->getX().at(0), p->getX().at(1), -bounds.at(2)};
                std::array<double, 3> v_arg = {p->getV().at(0), p->getV().at(1), -1*(p->getV().at(2))};
                double m_arg = p->getM();
                int type_arg = p->getType();
                Particle haloParticle = Particle(x_arg, v_arg, m_arg, type_arg);
                std::array<double, 3> force = calc.calculateLJF(&haloParticle, p, r_cutoff);
                std::array<double, 3> addedForce;
                for(int k = 0; k < 3; k++) {
                  addedForce[k] = p->getF()[k] + force[k];
                }
                p->setF(addedForce);
              } else { // i == 5      //Boundary to other XY Plane
                std::array<double, 3> x_arg = {p->getX().at(0),p->getX().at(1), cell_size.at(2)*cell_count.at(2) + bounds.at(i)};
                std::array<double, 3> v_arg = {p->getV().at(0), -p->getV().at(1), -1*(p->getV().at(2))};
                double m_arg = p->getM();
                int type_arg = p->getType();

                Particle haloParticle = Particle(x_arg, v_arg, m_arg, type_arg);
                std::array<double, 3> force = calc.calculateLJF(&haloParticle, p, r_cutoff);
                std::array<double, 3> addedForce;
                for(int k = 0; k < 3; k++) {
                  addedForce[k] = p->getF()[k] + force[k];
                }
                p->setF(addedForce);
              }
              //Particle has same absolute velocity but
              //eg. boundary is XY, then velocity in x and y are same but in z is -z;
              //Particle has same x and y position but z is same distance from boundary but reversed
              //Particle has no force

            } // if(bounds.at(i) != -1.0 && bounds.at(i) <= (pow(2, 1/6)/2)){

          } // for(int i = 0; i < 6; i++){

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
