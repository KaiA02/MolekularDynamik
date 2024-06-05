#ifndef CELL_H
#define CELL_H
#include <vector>
#include "../Particle.h"

class Cell {
private:
    std::vector<Particle*> particles; // Change to a vector of Particle pointers
    std::array<int, 3> id;
public:
    Cell(std::array<int, 3> id);
    std::vector<Particle*>& getParticles(); // Return a reference to the vector of Particle pointers
    void addParticle(Particle* p); // Change the parameter to a Particle pointer
    std::array<int, 3> getId();
    void emptyCell();
    bool isEmpty();
};
#endif //CELL_H
