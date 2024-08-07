//
// Created by jh on 29.05.2024.
//

#include "Cell.h"
#include "spdlog/spdlog.h"
Cell::Cell(){};
Cell::Cell(std::array<int, 3> init_id, bool halo) {
    id = init_id;
    halo_cell = halo;

}

bool Cell::isHalo() {
    return halo_cell;
}


std::vector<Particle*>& Cell::getParticles() {
    return particles;
}
std::array<int, 3> Cell::getId() {
    return id;
}


void Cell::addParticle(Particle* p) {
    particles.push_back(p);
}

void Cell::emptyCell() {
    particles.clear();
}

bool Cell::isEmpty() {
    return particles.empty();
}



