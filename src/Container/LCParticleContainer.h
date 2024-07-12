//
// Created by jh on 14.06.2024.
//

#ifndef LCPARTICLECONTAINER_H
#define LCPARTICLECONTAINER_H
#include "../Calculations.h"
#include "ParticleContainer.h"
#include "../utils/EpsilonSigma.h"
#include "../Membrane/Membrane.h"

/**
 * @brief The LCParticleContainer class
 * This class is a container for Particles using the linked cell algorithm.
 */
class LCParticleContainer : public ParticleContainer {
public:
 /**
 * @brief sets the cutoff
 */
  void setR_cutoff(double r_cutoff);

 /**
  * @brief sets the g_grav
  */
  void setG_grav(double g);

 /**
  * @brief getter for r_cutoff
  */
  double getR_cutoff();

 /**
  * @brief getter for g_grav
  */
  double getG_grav();

  /**
   *@brief realocates the particles to their new cells and also
   *implements the outflow boundary (managing the deltion of particles out of
   *the domain)
   */
  void realocateParticles();

  /**
   *@brief adds the particles to the cells
   */
  void fillCellsWithParticles();

  /**
   *@brief returns all particles of neighbouring cells
   *@param id id of the cells
   *@return a vector of particles which are in the neighbourhood of the cell
   *with id
   */
  std::vector<Particle> getParticleInNeighbourhood(std::array<int, 3> id);

  /**
   * @brief gets the cell with the according id
   * @param id of the cell
   * @return pointer to the cell
   * */
  Cell *getCellById(std::array<int, 3> id);
  /**
   *@brief generate all cells with according size
   *@param size_x size of x axis
   *@param size_y size of y axis
   *@param size_z size of z axis
   *@param r_cutoff  cutoff value for cells
   */
  void generateCells(int size_x, int size_y, int size_z, double r_cutoff);

  /**
   *@brief sets the boundary types(overflow, reflective, periodic, etc. )
   *@param in array of boundary types
   */
  void setBoundarys(std::array<int, 6> in);
  /**
   *@brief calls realocateParticle() and calculates LJF for all cells
   */
  void handleLJFCalculation(Calculations &calc, int timestep);
  /**
   *@brief adds the particle to all particles
   *and calls addParticleToCell(p)
   */
  void addParticle(Particle p);
  /**
   *@brief adds the particle to the according cell
   *returns true, if particle is correct added
   */
  bool addParticleToCell(Particle &p);

  /**
   * @return all particles
   */
  std::vector<Particle> &getParticles();

  /**
   *@brief adds all new Particles to the existing particle list
   *@param newParticles vector of new Particles
   */
  void addMultipleParticles(std::vector<Particle> &newParticles);

  /**
   * @return all cells
   */
  //std::vector<std::vector<std::vector<Cell>>> getCells();

  /**
   *@brief checks if the cell exists
   *@param id from according cell
   */
  bool cellExists(std::array<int, 3> id);

  /**
   *@brief only for programmers to debug
   */
  void countParticlesInCells();

  std::vector<Particle *> getBoundaryParticles();

  /**
   * @brief handles the boundary action according to the boundary type
   */
  void handleBoundaryAction();

  /**
   * @brief gets the distance of the particle to the boundaries
   * @param Particle the particle
   * @return an array with the distance to each boundary
   */
  std::array<double, 12> getInfluencingBoundarysWithDistance(Particle *);

  /**
   * @brief finds the position of the particle on the opposing side of the
   * domain (used for periodic boundary)
   * @param p the particle
   * @return the position of the particle on the opposing side of the domain
   */
  std::array<double, 3> findOponentXYZ(Particle *p);

  /**
   * finds the cell ID of the cell on the opposing side of the domain
   * @param ID the cell ID
   * @return the cell ID of the cell on the opposing side of the domain
   */
  std::array<int, 3> findOponentCellID(std::array<int, 3> ID);

  /**
   * calculates the force on the particle p using a halo particle
   * @param p the particle
   * @param x_arg the position of the halo particle
   */
  void calcWithHalo(Particle *p, std::array<double, 3> x_arg);

  /**
   * @brief applies the gravitational force to all particles
   */
  void applyGravitation();

  void setUpEpsilonAndSigmas();

  void setMembrane(Membrane m);

  Membrane getMembrane();

  std::vector<Particle*> getMovingParticles(std::array<std::array<int, 2>, 4> ids, int size);

 std::vector<Particle*> getAllParticlePointers();

private:
  std::vector<std::vector<std::vector<Cell>>> cells;
  std::array<double, 3> cell_size;
  std::array<int, 3> cell_count;
  std::array<int, 6> boundary_types;
  double r_cutoff;
  double g_grav;
  bool LJORSmoothLJ;
  std::vector<EpsilonSigma> epsAndSigs;
  Membrane membrane;

};

#endif // LCPARTICLECONTAINER_H
