#ifndef CELL_H
#define CELL_H
#include <vector>
#include "../Particle.h"

/**
   *@brief Cell class for generating a basic cell
   * used for the Linked-Cell algorithm in LCParticleContainer
   */
class Cell {
private:
    std::vector<Particle*> particles;
    std::array<int, 3> id;
    bool halo_cell;
    Cell* oposition;
public:
    Cell(std::array<int, 3> id, bool halo);
    void setOposition(Cell* op);
    Cell* getOposition();
    bool isHalo();
    std::vector<Particle*>& getParticles();
    /**
   *@brief adds particles to cell
   */
    void addParticle(Particle* p);
    std::array<int, 3> getId();
    /**
   *@brief clears all cells
   */
    void emptyCell();
    bool isEmpty();
};
#endif //CELL_H
