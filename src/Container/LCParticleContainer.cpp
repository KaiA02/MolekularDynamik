//
// Created by jh on 14.06.2024.
//

#include "LCParticleContainer.h"

#include "spdlog/spdlog.h"
#include "../Calculations.h"

void LCParticleContainer::setBoundarys(std::array<int, 6> in) {
  boundary_types = in;
}

double LCParticleContainer::getR_cutoff() {
  return r_cutoff;
}


void LCParticleContainer::setR_cutoff(double radius) {
  r_cutoff = radius;
}

double LCParticleContainer::getG_grav() {
  return g_grav;
}

void LCParticleContainer::setG_grav(double g) {
  g_grav = g;
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
  spdlog::debug("neighbourhood has {} particles", particleCount);
  if (particleCount <= 0) {
    return {};
  } else {
    spdlog::debug("there was a neigbourhood of {} particles", particleCount);
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
  spdlog::debug("we have {} cells and {} particles", cells.size(),getParticles().size());
  for (auto &c : cells) {
    c.emptyCell();
  }
  for (auto &p : particles) {
    int x = floor(p.getX().at(0) / cell_size.at(0));
    int y = floor(p.getX().at(1) / cell_size.at(1));
    int z = floor(p.getX().at(2) / cell_size.at(2));
    p.setF({0.0, 0.0, 0.0});
    if (cellExists({x, y, z})) {
      getCellById({x,y,z})->addParticle(&p);
    } else { //Parking of deleted Particles
        p.park();
      }
    }
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
      for (int x = -1; x <= cell_count[0]; x++) {
        for (int y = -1; y <= cell_count[1]; y++) {
          for (int z = -1; z <= cell_count[2]; z++) {
            bool isHalo;
            if(x == -1 || y== -1 || z == -1 || x == cell_count[0] || y == cell_count[1] || z == cell_count[2]) {
              isHalo = true;
              Cell c({x,y,z}, isHalo);
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
      spdlog::warn(
          "negative value detected in: generate Cells(count_x etc.) {}, {}, {}",
          cell_count[0], cell_count[1], cell_count[2]);
    }
  } else {
    spdlog::warn("negative value detected in: generate Cells(r_cutoff) {}",
                 r_cutoff);
  }
}

void LCParticleContainer::handleLJFCalculation(Calculations& calc) {
  realocateParticles();
  for (auto &c : cells) {
    if (!c.isEmpty()) {
      std::vector<Particle> neighbourhood = getParticleInNeighbourhood(c.getId());
      calc.LCcalculateLJF(c.getParticles(), neighbourhood);
    }
  }
  handleBoundaryAction();
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

void LCParticleContainer::countParticlesInCells() { //just for debugging
  int counter = 0;
  for(auto c: cells) {
    for(auto p: c.getParticles()) {
      counter ++;
    }
  }
  spdlog::info("there are {} particles in all the cells", counter);
}

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
      std::array<double, 12> bounds = getInfluencingBoundarysWithDistance(p);
        for(int i = 0; i < 6; i++) {
          if(boundary_types[i] == 1) {
            if(bounds.at(i) < 0) {
              p->park();
            }
          } //outflow
          else if(boundary_types[i] == 2){
            // generate halo Particle
            if(bounds.at(i) != -1.0 && abs(bounds.at(i)) <= pow(2, 1/6)){ //&& abs(bounds.at(i)) <= r_cutoff
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
            } // if(bounds.at(i) != -1.0 && bounds.at(i) <= (pow(2, 1/6)/2)){
          } //reflectiv
          else if(boundary_types[i] == 3) {
            if(bounds.at(i+6) >= 0.0) {
              std::array<double,3> x = findOponentXYZ(p);
              p->setX(x);
            }
          } //periodic

          }// if(boundary_types[i] == 2){

      }// for(auto p: boundaryparticles) {

  }

std::array<double, 12> LCParticleContainer::getInfluencingBoundarysWithDistance(Particle * p) {
  //returns distance to influencing Boundarys or -1 in case this Boundary isnt influencing. next 6 digits show, if Particle is ascending the Border(1) or descending (-1)
  std::array<double, 12> result = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  //get position and decide which boundary is affecting the Particle in the future
  double x = p->getX().at(0);
  double y = p->getX().at(1);
  double z = p->getX().at(2);

  int x_id = int(floor(p->getX().at(0) / cell_size.at(0)));
  int y_id = int(floor(p->getX().at(1) / cell_size.at(1)));
  int z_id = int(floor(p->getX().at(2) / cell_size.at(2)));

  if(x_id == 0){ //Particle close to boundary 1
    double distance1 = x;
    result.at(0) = distance1;
    if(p->getV().at(0) < 0.0) {
      result.at(6) = 1.0;
    } else if(p->getV().at(0) > 0.0) {
      result.at(6) = -1.0;
    }
  }
  if(x_id == cell_count.at(0) -1) { //Particle close to boundary 2
    double distance2 = cell_size.at(0) - (x - (cell_size.at(0)*(cell_count.at(0)-1)));
    result.at(1) = distance2;
    if(p->getV().at(0) < 0.0) {
      result.at(7) = -1.0;
    } else if(p->getV().at(0) > 0.0) {
      result.at(7) = 1.0;
    }
  }
  if(y_id == 0){ //Particle close to boundary 3
    double distance3 = y;
    result.at(2) = distance3;
    if(p->getV().at(1) < 0.0) {
      result.at(8) = 1.0;
    } else if(p->getV().at(1) > 0.0) {
      result.at(8) = -1.0;
    }
  }
  if(y_id == cell_count.at(1) -1) { //Particle close to boundary 4
    double distance4 = cell_size.at(1) - (y - (cell_size.at(1)*(cell_count.at(1)-1)));
    result.at(3) = distance4;
    if(p->getV().at(1) > 0.0) {
      result.at(9) = 1.0;
    } else if(p->getV().at(1) < 0.0) {
      result.at(9) = -1.0;
    }
  }
  if(z_id == 0){ //Particle close to boundary 5
    double distance5 = z;
    result.at(4) = distance5;
    if(p->getV().at(2) < 0.0) {
      result.at(10) = 1.0;
    } else if(p->getV().at(2) > 0.0) {
      result.at(10) = -1.0;
    }
  }
  if(z_id >= cell_count.at(2) -1) { //Particle close to boundary 6
    double distance6 = cell_size.at(2) - (z - (cell_size.at(2)*(cell_count.at(2)-1)));
    result.at(5) = distance6;
    if(p->getV().at(2) > 0.0) {
      result.at(11) = 1.0;
    } else if(p->getV().at(2) < 0.0) {
      result.at(11) = -1.0;
    }
  }

  return result;
}

std::array<double, 3> LCParticleContainer::findOponentXYZ(Particle* p) {

  //calculates Cell ID
  int x_id = floor(p->getX().at(0) / cell_size.at(0));
  int y_id = floor(p->getX().at(1) / cell_size.at(1));
  int z_id = floor(p->getX().at(2) / cell_size.at(2));

  // finds oponentCell ID
  std::array<int, 3> oponentCellID = findOponentCellID({x_id, y_id, z_id});

  double x_cell = x_id * cell_size.at(0);
  double y_cell = y_id * cell_size.at(1);
  double z_cell = z_id * cell_size.at(2);

  double x_relative_to_cell = p->getX().at(0) - x_cell;
  double y_relative_to_cell = p->getX().at(1) - y_cell;
  double z_relative_to_cell = p->getX().at(2) - z_cell;

  //calcs right XYZ for Particle in new Cell
  double x_oponentCell = oponentCellID.at(0) * cell_size.at(0);
  double y_oponentCell = oponentCellID.at(1) * cell_size.at(1);
  double z_oponentCell = oponentCellID.at(2) * cell_size.at(2);

  double x_relative_to_OponentCell = x_oponentCell + x_relative_to_cell;
  double y_relative_to_OponentCell = y_oponentCell + y_relative_to_cell;
  double z_relative_to_OponentCell = z_oponentCell + z_relative_to_cell;

  //returns XYZ of Particle in new Cell
  return {x_relative_to_OponentCell, y_relative_to_OponentCell, z_relative_to_OponentCell};
}

std::array<int, 3> LCParticleContainer::findOponentCellID(std::array<int, 3> ID) {
  std::array<int, 3> return_id = ID;
  if(ID.at(0) < 0) {
    return_id.at(0) = cell_count.at(0) -1;
  } else if(ID.at(0) >= cell_count.at(0)) {
    return_id.at(0) = 1;
  }
  if(ID.at(1) < 0) {
    return_id.at(1) = cell_count.at(1) -1;
  } else if(ID.at(1) >= cell_count.at(1)) {
    return_id.at(1) = 1;
  }
  if(ID.at(2) < 0) {
    return_id.at(2) = cell_count.at(2) -1;
  } else if(ID.at(2) >= cell_count.at(2)) {
    return_id.at(2) = 1;
  }
  return return_id;
}

void LCParticleContainer::calcWithHalo(Particle *p, std::array<double, 3> x_arg, std::array<double, 3> v_arg) {
  double m_arg = p->getM();
  int type_arg = p->getType();
  Particle haloParticle = Particle(x_arg, v_arg, m_arg, type_arg);
  Calculations calc(*this);
  std::array<double, 3> addedForce;
  std::array<double, 3> f_ij = calc.calculateLJF(p, &haloParticle);
  for(int k = 0; k < 3; k++) {
    addedForce[k] = p->getF()[k] + f_ij[k];
  }
  p->setF(addedForce);
}

void LCParticleContainer::applyGravitation() {
  for(auto p : particles) {
    double y_force = p.getX().at(1);
    y_force += g_grav * p.getM();
    p.setX({p.getX().at(0), y_force, p.getX().at(2)});
  }
}
