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
    ParticleGenerator();
    void generateCuboid(const Particle& start, int n1, int n2, int n3, double distance, double meanVelocity);
    std::vector<Particle> getCube();

};

#endif //PARTICLEGENERATOR_H
