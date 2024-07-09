//
// Created by jh on 04.07.2024.
//

#ifndef MEMBRANEPAIR_H
#define MEMBRANEPAIR_H

#include "../Particle.h"
#include <vector>

class MembranePair{
    public:
        MembranePair(Particle* center);
        Particle* getCenter();
        void setCenter(Particle* center);
        std::vector<Particle*> getDirect();
        void appendDirect(Particle* p);
        std::vector<Particle*> getDiagonal();
        void appendDiagonal(Particle* p);
    private:
        Particle* center;
        std::vector<Particle*> direct;
        std::vector<Particle*> diagonal;
};
#endif //MEMBRANEPAIR_H
