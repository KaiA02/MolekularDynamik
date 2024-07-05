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


std::vector<Particle> LCParticleContainer::getParticleInNeighbourhood(std::array<int, 3> id) {
  std::vector<Particle> neigbourhood;
  int particleCount = 0;
  std::array<int, 3> possible_x = {id[0]-1, id[0], id[0]+1};
  std::array<int, 3> possible_y = {id[1]-1, id[1], id[1]+1};
  std::array<int, 3> possible_z = {id[2]-1, id[2], id[2]+1};
  Cell* cell = nullptr;
  for(int x = 0; x < 3; x++){
    for(int y = 0; y < 3; y++){
      for(int z = 0; z < 3; z++){
        if(x!=1 || y!=1 || z!=1){
          if(cellExists({possible_x[x], possible_y[y], possible_z[z]})){
            cell = getCellById({possible_x[x], possible_y[y], possible_z[z]});
            for(auto p: cell->getParticles()){
              neigbourhood.push_back(*p);
              particleCount++;
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

Cell* LCParticleContainer::getCellById(const std::array<int, 3> id) {
  spdlog::debug("getCellByID: id is {} {} {}", id[0]+1, id[1]+1, id[2]+1);
  return &cells.at(id[0]+1).at(id[1]+1).at(id[2]+1);
}

void LCParticleContainer::realocateParticles() {
  std::vector<std::vector<std::vector<Cell>>> newCells = std::vector<std::vector<std::vector<Cell>>>{};
  int x;
  int y;
  int z;
  // cells size in x is set to cell_count in x plus 2(halo_cells)
  newCells.resize(cell_count[0]+2);

  // cells size in y is set to cell_count in y plus 2(halo_cells)
  for (int i = 0; i < cell_count[0]+2; ++i) {
    newCells[i].resize(cell_count[1]+2);
  }

  // cells size in z is set to cell_count in z plus 2(halo_cells)
  for (int i = 0; i < cell_count[0]+2; ++i) {
    for (int j = 0; j < cell_count[1]+2; ++j) {
      newCells[i][j].resize(cell_count[2]+2);
    }
  }
  cells = newCells;
  for (auto &p : particles) {
    spdlog::debug("RelaocateParticle: will realocate Particle");
    x = floor(p.getX().at(0) / cell_size.at(0));
    y = floor(p.getX().at(1) / cell_size.at(1));
    z = floor(p.getX().at(2) / cell_size.at(2));
    p.setF({0.0, 0.0, 0.0});
    if (cellExists({x, y, z})) {
      spdlog::debug("RelaocateParticle: cell exists");
      getCellById({x,y,z})->addParticle(&p);
      spdlog::debug("RelaocateParticle: added Particle to Cell");
    } else { //Parking of deleted Particles
      //if(x > -20 && y > -20 && z > -20){
      //  spdlog::debug("RelaocateParticle: Cell id {} {} {} dosent exists", x, y, z);
      //}

      p.park();

      }
    }
    spdlog::debug("all particles are added--------------------------------------->");
  }

void LCParticleContainer::fillCellsWithParticles() {
  int counter =0;
  int addedCounter = 0;
  for (auto &p : particles) {
    bool added = addParticleToCell(p);
    counter ++;
    if(added) {
      addedCounter++;
    }
  }
  spdlog::info("called addParticleToCell() {} times and it was {} times succesfull", counter, addedCounter);
}

void LCParticleContainer::generateCells(const int size_x, const int size_y, const int size_z, const double r_cutoff) {
  std::vector<std::vector<std::vector<Cell>>> newCells{};
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

      // cells size in x is set to cell_count in x plus 2(halo_cells)
      newCells.resize(cell_count[0]+2);

      // cells size in y is set to cell_count in y plus 2(halo_cells)
      for (int i = 0; i < cell_count[0]+2; ++i) {
        newCells[i].resize(cell_count[1]+2);
      }

      // cells size in z is set to cell_count in z plus 2(halo_cells)
      for (int i = 0; i < cell_count[0]+2; ++i) {
        for (int j = 0; j < cell_count[1]+2; ++j) {
          newCells[i][j].resize(cell_count[2]+2);
        }
      }

      spdlog::info("Cellcounts: {} {} {}", cell_count[0], cell_count[1], cell_count[2]);
      spdlog::info("Cellsize: {} {} {}", cell_size_x, cell_size_y, cell_size_z);
      spdlog::info("Domainsize: {} {} {}", size_x, size_y, size_z);
      spdlog::info("r_cutoff: {}", r_cutoff);
      // Clear existing cells before generating new ones
      for (int x = -1; x <= cell_count[0]; x++) {
        for (int y = -1; y <= cell_count[1]; y++) {
          for (int z = -1; z <= cell_count[2]; z++) {
            bool isHalo = false;
            if(x == -1 || y== -1 || z == -1 || x == cell_count[0] || y == cell_count[1] || z == cell_count[2]) {
              isHalo = true;
            }
            Cell c({x,y,z}, isHalo);
            newCells.at(x+1).at(y+1).at(z+1) = c;
          }
        }
      }
      cell_size = {cell_size_x, cell_size_y, cell_size_z};
      cells = newCells;
      spdlog::debug("generated Cells");
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

void LCParticleContainer::handleLJFCalculation(Calculations& calc, int timestep) {
  realocateParticles();
  spdlog::debug("LCPartCon: realoctaed Particles");
  std::vector<Particle> neighbourhood;
  Cell* c = nullptr;
  spdlog::debug("LCPartCon: itterating over all neighbourhoods");
  for(std::vector<Cell>::size_type x = 0; x < cells.size() -1; x++){
    for(std::vector<Cell>::size_type y = 0; y < cells[x].size() -1; y++){
      for(std::vector<Cell>::size_type z = 0; z < cells[x][z].size() -1; z++){
         c = getCellById({int(x), int(y), int(z)});
          if (!c->isEmpty()) {
            neighbourhood = getParticleInNeighbourhood({int(x), int(y), int(z)});
            calc.LCcalculateLJF(c->getParticles(), neighbourhood, epsAndSigs);
          }
        }
      }
    }
  handleBoundaryAction();
  applyGravitation();
  if(timestep <= 150){
    membrane.applyMovement();
  }
  membrane.stabilizeMembrane(calc);

}

bool LCParticleContainer::addParticleToCell(Particle &p) {
  const int x = int(floor(p.getX().at(0) / cell_size.at(0)));
  const int y = int(floor(p.getX().at(1) / cell_size.at(1)));
  const int z = int(floor(p.getX().at(2) / cell_size.at(2)));
  if (cellExists({x, y, z})) {
    getCellById({x, y, z})->addParticle(&p);
    return true;
  }
  return false;
}

std::vector<Particle> &LCParticleContainer::getParticles() { return particles; }

void LCParticleContainer::addMultipleParticles(std::vector<Particle> &newParticles) {
  int counter = 0;
  for (auto &p : newParticles) {
    particles.push_back(p);
    counter++;
  }
  spdlog::info("added {} particles from generator to Container", counter);
}

void LCParticleContainer::addParticle(Particle p) {
  particles.push_back(p);
  addParticleToCell(p);
}

//std::vector<std::vector<std::vector<Cell>>> LCParticleContainer::getCells() { return cells; }

bool LCParticleContainer::cellExists(std::array<int, 3> id) {
  spdlog::debug("cell exists: id is {} {} {}", id[0], id[1], id[2]);
  return(id[0] >= -1 && id[0] <= cell_count[0]  && id[1] >= -1 && id[1] <= cell_count[1]  && id[2] >= -1 && id[2] <= cell_count[2]);
};



std::vector<Particle *> LCParticleContainer::getBoundaryParticles() {
  std::vector<Particle*> result;
  const std::vector<int> possible_x = {-1, 0, cell_count[0] -1 , cell_count[0]};
  const std::vector<int> possible_y = {-1, 0, cell_count[1] -1 , cell_count[1]};
  std::vector<int> possible_z;
  if(cell_count.at(2) >= 2) {
    possible_z = {-1, 0, cell_count[2] -1 , cell_count[2]};
  } else {
    possible_z = {-1, 0, cell_count[2]};
  }
  //runtime is now 4 * CellCount[1] * CellCount[2] + CellCount[0] * 4 * CellCount[2] + CellCount[0] * CellCount[1] * 4
  // instead of CellCount[0] * CellCount[1] * CellCount[2] e.g CellCount[0]=16 CellCount[1]=16 CellCount[2]=16
  //then it was 4096 an now it is 3072. So with big Domains and higher  CellCounts in mind this is faster
  for(auto x : possible_x) {
    for(int y = -1; y < cell_count[1] +1; y++) {
      for(int z = -1; z < cell_count[2] +1; z++) {
        if (cellExists({x,y,z})) {
          for(auto p : getCellById({x,y,z})->getParticles()) {
            result.push_back(p);
          }
        }
      }
    }
  }
  for(int x = -1; x < cell_count[0] +1; x++) {
    for(auto y : possible_y) {
      for(int z = -1; z < cell_count[2] +1; z++) {
        if (cellExists({x,y,z})) {
          for(auto p : getCellById({x,y,z})->getParticles()) {
            result.push_back(p);
          }
        }
      }
    }
  }
  for(int x = -1; x < cell_count[0] +1; x++) {
    for(int y = -1; y < cell_count[1] +1; y++) {
      for(auto z : possible_z) {
        if (cellExists({x,y,z})) {
          for(auto p : getCellById({x,y,z})->getParticles()) {
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
  std::array<double, 12> bounds{};
  std::array<double, 3> X{};
    for(auto p: boundaryparticles) {
      X = p->getX();
      bounds = getInfluencingBoundarysWithDistance(p);
        for(int i = 0; i < 6; i++) {
          if(boundary_types[i] == 1) {
            if(bounds.at(i) < 0) {
              p->park();
            }
          } //outflow
          else if(boundary_types[i] == 2){
            // generate halo Particle
            if(bounds.at(i) != -1.0 && abs(bounds.at(i)) <= 1.122462048309373){ //&& abs(bounds.at(i)) <= r_cutoff
              if( i == 0 ){ //Boundary to YZ Plane
                std::array<double, 3> x_arg = {-bounds.at(i), X.at(1), X.at(2)};
                calcWithHalo(p, x_arg);
              } else if ( i == 1 ) { //Boundary to other YZ Plane
                std::array<double, 3> x_arg = {cell_size.at(0)*cell_count.at(0) + bounds.at(i), X.at(1), X.at(2)};
                calcWithHalo(p, x_arg);
              } else if ( i == 2 ) { //Boundary to XZ Plane
                std::array<double, 3> x_arg = {X.at(0), -bounds.at(i), X.at(2)};
                calcWithHalo(p, x_arg);
              } else if ( i == 3 ) { //Boundary to other XZ Plane
                std::array<double, 3> x_arg = {X.at(0), cell_size.at(1)*cell_count.at(1) + bounds.at(i), X.at(2)};
                calcWithHalo(p, x_arg);
              } else if ( i == 4 ) { //Boundary to XY Plane
                std::array<double, 3> x_arg = {X.at(0), X.at(1), -bounds.at(i)};
                calcWithHalo(p, x_arg);
              } else { // i == 5      //Boundary to other XY Plane
                std::array<double, 3> x_arg = {X.at(0),X.at(1), cell_size.at(2)*cell_count.at(2) + bounds.at(i)};
                calcWithHalo(p, x_arg);
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
  const std::array<double, 3> xyz = p->getX();
  const double x = xyz.at(0);
  const double y = xyz.at(1);
  const double z = xyz.at(2);

  const std::array<double ,3> V = p->getV();

  const int x_id = int(floor(x / cell_size.at(0)));
  const int y_id = int(floor(y / cell_size.at(1)));
  const int z_id = int(floor(z / cell_size.at(2)));

  if(x_id == 0){ //Particle close to boundary 1
    result.at(0) = x;
    if(V[0] < 0.0) {
      result.at(6) = 1.0;
    }
  }
  if(x_id == cell_count.at(0) -1) { //Particle close to boundary 2
    result.at(1) = cell_size.at(0) - (x - (cell_size.at(0)*(cell_count.at(0)-1)));
    if(V[0] > 0.0) {
      result.at(7) = 1.0;
    }
  }
  if(y_id == 0){ //Particle close to boundary 3
    result.at(2) = y;
    if(V[1] < 0.0) {
      result.at(8) = 1.0;
    }
  }
  if(y_id == cell_count.at(1) -1) { //Particle close to boundary 4
    result.at(3) = cell_size.at(1) - (y - (cell_size.at(1)*(cell_count.at(1)-1)));
    if(V[1] > 0.0) {
      result.at(9) = 1.0;
    }
  }
  if(z_id == 0){ //Particle close to boundary 5
    result.at(4) = z;
    if(V[2] < 0.0) {
      result.at(10) = 1.0;
    }
  }
  if(z_id >= cell_count.at(2) -1) { //Particle close to boundary 6
    result.at(5) = cell_size.at(2) - (z - (cell_size.at(2)*(cell_count.at(2)-1)));
    if(V[2] > 0.0) {
      result.at(11) = 1.0;
    }
  }
  return result;
}

std::array<double, 3> LCParticleContainer::findOponentXYZ(Particle* p) {

  //calculates Cell ID
  const int x_id = floor(p->getX().at(0) / cell_size.at(0));
  const int y_id = floor(p->getX().at(1) / cell_size.at(1));
  const int z_id = floor(p->getX().at(2) / cell_size.at(2));

  // finds oponentCell ID
  std::array<int, 3> oponentCellID = findOponentCellID({x_id, y_id, z_id});

  //calcs the rest
  return {oponentCellID.at(0) * cell_size.at(0) + p->getX().at(0) - x_id * cell_size.at(0),
            oponentCellID.at(1) * cell_size.at(1) + p->getX().at(1) - y_id * cell_size.at(1),
            oponentCellID.at(2) * cell_size.at(2) + p->getX().at(2) - z_id * cell_size.at(2)};
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

void LCParticleContainer::calcWithHalo(Particle *p, std::array<double, 3> x_arg) {
  Particle haloParticle = Particle(x_arg, {}, p->getM(), p->getType());
  Calculations calc(*this);
  std::array<double, 3> addedForce;
  std::array<double, 3> f_ij = calc.calculateLJF(p, &haloParticle, p->getEpsilon(), p->getSigma());
  addedForce = {p->getF().at(0) + f_ij.at(0), p->getF().at(1) + f_ij.at(1), p->getF().at(2) + f_ij.at(2)};
  p->setF(addedForce);
}

void LCParticleContainer::applyGravitation() {
  std::array<double, 3> F;
  for(auto &p : particles) {
    F = p.getF();
    p.setF({F[0], F[1] + (g_grav * p.getM()), F[2]});
  }
}

void LCParticleContainer::setUpEpsilonAndSigmas() {
  int type = 0;
  std::vector<EpsilonSigma> oldOne {};
  bool found;
  for(auto p : particles) {
    type = p.getType();
    found = false;
    for(auto es : oldOne) {
      if(!found) {
        if(es.getT1() == type) {
          found = true;
        }
      }
    }
    if(!found) {
      EpsilonSigma PWithP (type, type, p.getSigma(), p.getEpsilon());
    oldOne.push_back(PWithP);
    }
  }
  std::vector<EpsilonSigma> newOne {};
  int i_type;
  double i_epsilon;
  double i_sigma;

  for(std::vector<EpsilonSigma>::size_type i= 0; i < oldOne.size(); i++){

    i_epsilon = oldOne.at(i).getEpsilon();
    i_sigma = oldOne.at(i).getSigma();
    i_type = oldOne.at(i).getT1();

    for(std::vector<EpsilonSigma>::size_type j= 0; j < oldOne.size(); j++){
      double sigma = (i_sigma +  oldOne.at(j).getSigma()) * 0.5;
      double epsilon = std::sqrt(i_epsilon * oldOne.at(j).getEpsilon());
      EpsilonSigma IWithJ (i_type, oldOne.at(j).getT1(), sigma, epsilon);
      newOne.push_back(IWithJ);
    }
    epsAndSigs = newOne;
  }
}

void LCParticleContainer::setMembrane(Membrane m){
  membrane = m;
}
Membrane LCParticleContainer::getMembrane(){
  return membrane;
};

