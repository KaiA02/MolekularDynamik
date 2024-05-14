//
// Created by jh on 13.05.2024.
//

#ifndef PARTICLEGENERATOR_H
#define PARTICLEGENERATOR_H
#include <vector>

#include "Particle.h"

class ParticleGenerator {
private:
    std::vector<Particle> cube;

public:
    ParticleGenerator(const Particle& start, int sizeX, int sizeY, int sizeZ, double distance, double meanVelocity);
    std::vector<Particle> getCube();

};

#endif //PARTICLEGENERATOR_H
