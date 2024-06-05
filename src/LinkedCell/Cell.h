//
// Created by jh on 29.05.2024.
//

#ifndef CELL_H
#define CELL_H
#include <vector>
#include "../Particle.h"

class Cell {
private:
    std::vector<Particle> particles;
    std::array<int, 3> id;
public:
    Cell(std::array<int, 3> id);
    std::vector<Particle>& getParticles();
    void addParticle(const Particle& p);
    std::array<int, 3> getId();
    void emptyCell();
    bool isEmpty();
};
#endif //CELL_H
