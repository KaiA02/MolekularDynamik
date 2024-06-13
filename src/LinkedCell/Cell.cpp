//
// Created by jh on 29.05.2024.
//

#include "Cell.h"
#include "spdlog/spdlog.h"

Cell::Cell(std::array<int, 3> init_id) {
    id = init_id;
}

std::vector<Particle*>& Cell::getParticles() {
    return particles;
}
std::array<int, 3> Cell::getId() {
    return id;
}


void Cell::addParticle(Particle* p) {
    size_t pre = particles.size();
    particles.push_back(p);
    size_t post = particles.size();
    //only for debug
    if(pre == post) {
        spdlog::debug("error ist the same when added--------------------------------");
    }
}

void Cell::emptyCell() {
    particles.clear();
}

bool Cell::isEmpty() {
    return particles.empty();
}



