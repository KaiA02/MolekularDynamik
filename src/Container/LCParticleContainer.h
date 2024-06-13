//
// Created by kaiarenja on 13.06.24.
//

#ifndef LCPARTICLECONTAINER_H
#define LCPARTICLECONTAINER_H
#include "ParticleContainer.h"


/**
 * @brief The LCParticleContainer class
 * This class is a container for Particles using the linked cell algorithm.
 */
class LCParticleContainer : public ParticleContainer {
public:
 LCParticleContainer();
 /**
  *@brief realocates the particles to their new cells and also
  *implements the outflow boundary (managing the deltion of particles out of the domain)
 */
 void realocateParticles();

 /**
  *@brief adds the particles to the cells
 */
 void fillCellsWithParticles();

 /**
  *@brief returns all particles of neighbouring cells
  *@param id id from the celss
 */
 std::vector<Particle*> getParticleInNeighbourhood(std::array<int, 3> id);

 Cell* getCellById(std::array<int, 3> id);
 /**
  *@brief generate all cells with according size
  *@param size_x size of x axis
  *@param size_y size of y axis
  *@param size_z size of z axis
  *@param r_cutoff  cutoff value for cells
 */
 void generateCells(int size_x, int size_y, int size_z, double r_cutoff);

 void setBoundarys(std::array<int, 6> in);
 /**
  *@brief calls realocateParticle() and calculates LJF for all cells
 */
 void handleLJFCalculation();
 /**
  *@brief adds the particle to all particles
  *and calls addParticleToCell(p)
 */
 void addParticle(Particle p);
 /**
 *@brief adds the particle to the according cell
 *returns true, if particle is correct added
*/
 bool addParticleToCell(Particle& p);
 std::vector<Particle>& getParticles();
 /**
 *@brief adds all new Particles to the existing particle list
 *@param newParticles vector of new Particles
*/
 void addMultipleParticles(std::vector<Particle>& newParticles);
 std::vector<Cell> getCells();
 /**
  *@brief checks if the cell exists
  *@param id from according cell
 */
 bool cellExists(std::array<int, 3> id);
 /**
  *@brief only for programmers to debug
 */
 void countParticlesInCells();


 std::vector<Particle*> getBoundaryParticles();
 void handleBoundaryAction();
 std::array<double, 6> getInfluencingBoundarysWithDistance(Particle*);

private:
 std::vector<Cell> cells;
 std::array<double, 3> cell_size;
 std::array<int, 3> cell_count;
 std::array<int, 6> boundary_types;


};

#endif //LCPARTICLECONTAINER_H