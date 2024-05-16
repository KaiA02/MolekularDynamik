//
// Created by kaiarenja on 16.05.24.
//

#ifndef CALCULATIONS_H
#define CALCULATIONS_H
#include "ParticleContainer.h"

class Calculations {
private:
    ParticleContainer& particles;
public:
    Calculations(ParticleContainer &other);
    /**
    * calculate the force for all particles
    */
    void calculateF();
    void calculateLJF();
    /**
    * calculate the position for all particles
    */
    void calculateX(double delta_t);
    /**
    * calculate the velocity for all particles
    */
    void calculateV(double delta_t);

    ParticleContainer& getParticles();
};

#endif //CALCULATIONS_H
