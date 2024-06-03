
#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include "Particle.h"
#include <vector>

#include "LinkedCell/Cell.h"
/**
 * @brief The ParticleContainer class
 * This class is a container for Particles. It contains a vector of Particles
 * and a 2D vector of pairings.
 */
class ParticleContainer {
public:
  ParticleContainer();
  /**
   * @brief add a Particle to the vector
   * @param Particle particle
   */
  void addParticle(const Particle &particle);

  /**
  @brief add a Pairing to the 2d-vector pairings
  @param index1: index of the first Particle in particles index2: index of the
  second Particle in particles
   */
  void addPairing(int particleIndex1, int particleIndex2);
  std::vector<Particle> &getParticles();
  const std::vector<Particle> &getParticles() const;
  int size() const;

  std::vector<Particle>::iterator begin();
  std::vector<Particle>::iterator end();
  std::vector<Particle>::const_iterator begin() const;
  std::vector<Particle>::const_iterator end() const;

  /**
   * @brief adds a cuboid of particles to the vector particles
   * @param particleCube cuboid generated by the ParticleGenerator
   */
  void addCube(std::vector<Particle> particleCube);

 /**
   * @brief adds a disk of particles to the vector particles
   * @param particleDisk disk generated by the ParticleGenerator
   */
 void addDisk(std::vector<Particle> particleDisk);

  /**
   * @brief replaces a particle in the vector particles at index position
   * @param p new Particle
   * @param position index in particles
   */
  void setParticle(Particle p, int position);

  void resetParticles();

private:
  std::vector<Particle> particles;

  // TODO: integrate pairings into calculations when number of particles gets
  // too big (maybe use a force threshold to determine which particles are close
  // enough to interact with each other)
  std::vector<std::pair<int, int>> pairings;
};
class LCParticleContainer: public ParticleContainer {
private:
 std::vector<Cell> cells;
 std::array<double, 3> cell_size;

public:
 //input CellSize in x,y,z and Strategy for handling particles out of border
 //allocates Particle to right cell and resets the forces to zero
 void realocateParticles(int handle_out_of_border);
 std::vector<Particle> getParticleInNeighbourhood(Cell c);
 Cell* getCellById(std::array<int, 3> id);
 //input: Domain Size in x,y,z and r_cutoff
 void generateCells(int size_x, int size_y, int size_z, double r_cutoff);
 //handles LJFCalcualtion for all Cells;
 void handleLJFCalculation();
 void addParticle(Particle p);
};

#endif // PARTICLECONTAINER_H
