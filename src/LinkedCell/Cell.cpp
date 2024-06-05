//
// Created by jh on 29.05.2024.
//

#include "Cell.h"
Cell::Cell(std::array<int, 3> init_id) {
    id = init_id;
}
std::array<int, 3> Cell::getId() {
    return id;
}
std::vector<Particle>& Cell::getParticles() {
    return particles;
}
void Cell::addParticle(const Particle& p) {
    particles.push_back(p);
}
void Cell::emptyCell() {
    particles = {};
}
bool Cell::isEmpty() {
    return particles.size() == 0;
}



