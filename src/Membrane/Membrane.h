//
// Created by jh on 04.07.2024.
//

#ifndef MEMBRANE_H
#define MEMBRANE_H
#include <vector>
#include "../Particle.h"
#include "MembranePair.h"
#include "../Calculations.h"
class Membrane{
    public:
        Membrane();
        Membrane(std::vector<Particle*> particles, double distance, double fUp);
        void applyMovement();
        void stabilizeMembrane(Calculations& calc);
        void setMovingParticles(std::vector<Particle*> movPart);
        int getMovingParticleCount();
        void getAveragePairSize();

    private:
        std::vector<MembranePair> pairs;
        double distance;
        std::vector<Particle*> movingParticles;
        double forceUpwards;
};
#endif //MEMBRANE_H
