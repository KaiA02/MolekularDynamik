#ifndef CELL_H
#define CELL_H
#include "../Particle.h"
#include <vector>

/**
 *@brief Cell class for generating a basic cell
 * used for the Linked-Cell algorithm in LCParticleContainer
 */
class Cell {
private:
  std::vector<Particle *> particles;
  std::array<int, 3> id;
  bool halo_cell;

public:
  Cell(std::array<int, 3> id, bool halo);

  /**
   * @brief sets the opposite cell of the current cell
   * @param op the opposite cell
   */
  void setOposition(Cell *op);

  /**
   *
   * @return the opposite cell
   */
  Cell *getOposition();

  /**
   *
   * @return true if the cell is a halo cell, false otherwise
   */
  bool isHalo();

  /**
   *
   * @return the particles in the cell
   */
  std::vector<Particle *> &getParticles();
  /**
   *@brief adds particles to cell
   */
  void addParticle(Particle *p);
  std::array<int, 3> getId();
  /**
   *@brief clears all cells
   */
  void emptyCell();

  /**
   *
   * @return true if the cell is empty, false otherwise
   */
  bool isEmpty();
};
#endif // CELL_H
