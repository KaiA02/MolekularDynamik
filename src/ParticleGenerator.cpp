//
// Created by jh on 13.05.2024.
//

#include "ParticleGenerator.h"
#include <utils/MaxwellBoltzmannDistribution.h>

#include <vector>

ParticleGenerator::ParticleGenerator(Particle start, int n1, int n2, int n3, double distance, double meanVelocity) {
    std::array<double, 3> maxwellVelocity = maxwellBoltzmannDistributedVelocity(meanVelocity, 3);
    for(int x = 0; x < n1; ++x) {
        for(int y = 0; y < n2; ++y) {
            for(int z = 0; z < n3; ++z) {
                Particle p(0);
                std::array<double, 3> newPosition = {start.getX()[0] + x * distance, start.getX()[1] + y * distance, start.getX()[2] + z * distance};
                p.setX(newPosition);
                std::array<double, 3> addedVelocity;
                for(int i = 0; i < 3; i++) {
                    addedVelocity[i] = start.getV()[i] + maxwellVelocity[i];
                }
                p.setV(addedVelocity);

                p.setF(start.getF());
                p.setM(start.getM());
                p.setType(start.getType());

                cube[x][y][z] = p;


            }
        }
    }
}

std::vector<std::vector<std::vector<Particle>>> ParticleGenerator::getCube() {
    return cube;
}
