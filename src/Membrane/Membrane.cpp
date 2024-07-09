//
// Created by jh on 04.07.2024.
//
#include "Membrane.h"
//#include "../Container/LCParticleContainer.h"
#include <cmath>
#include "spdlog/spdlog.h"

Membrane::Membrane(){
 distance = 0;
 std::vector<MembranePair> pair;
 pairs = pair;
 std::vector<Particle*> movingParticle;
 movingParticles = movingParticle;
forceUpwards= 0;
};

void Membrane::getAveragePairSize(){
int counter = 1;
    for(auto pair: pairs){
        spdlog::warn("pair id is {}:", counter);
        spdlog::warn("    has {} direct and {} diagonal", pair.getDirect().size(), pair.getDiagonal().size());
        spdlog::warn("-----------------------------");
        counter ++;
    }

}


Membrane::Membrane(std::vector<Particle*> particles, double d, double fUp){
    distance = d;
    forceUpwards = fUp;
    double dist;
    std::array<double,3> x1;
    std::array<double,3> x2;
    int directCounter = 0;
    int diagonalCounter = 0;
    int centerCounter = 0;
    for(auto p1:particles){
        centerCounter ++;
        x1 = p1->getX();
        MembranePair pair(p1);
        pair.setCenter(p1);
        for(auto p2:particles){
            x2 = p2->getX();
            dist = sqrt(pow((x1[0] - x2[0]),2) + pow((x1[1] - x2[1]),2) + pow((x1[2] - x2[2]),2));
            if(dist != 0 ){
                if(dist <= distance+0.001 && distance >= -distance-0.001){ // direct Neighbour
                    pair.appendDirect(p2);
                    directCounter ++;
                } else if(dist <= sqrt(distance*distance*2)+0.001 && dist >= -1*sqrt(distance*distance*2)-0.001){ //diagonal Neighbour
                    pair.appendDiagonal(p2);
                    diagonalCounter ++;
                }
            }
        }
        pairs.push_back(pair);
    }
    getAveragePairSize();
};

void Membrane::applyMovement(){
    std::array<double, 3> force;
    std::array<double, 3> newForce;
    for(auto p: movingParticles){
        force = p->getF();
        newForce = {force[0], force[1], force[2]+forceUpwards};
        p->setF(newForce);
    }
};

void Membrane::stabilizeMembrane(Calculations& calc){
    std::array<double,3> f_ij;
    Particle* center;
    std::array<double,3> centerForce;
    std::array<double,3> otherForce;
    int pair_count = 0;
    for(auto pair: pairs){
        center = pair.getCenter();
        centerForce = center->getF();
        for(auto dir: pair.getDirect()){
            otherForce = dir->getF();
            f_ij = calc.calculateHarmonicForce(center, dir);
            center->setF({centerForce[0]+f_ij[0], centerForce[1]+f_ij[1], centerForce[2]+f_ij[2]});
            dir->setF({otherForce[0]-f_ij[0], otherForce[1]-f_ij[1], otherForce[2]-f_ij[2]});
        }
        for(auto dia: pair.getDiagonal()){
            otherForce = dia->getF();
            f_ij = calc.calculateHarmonicForceDiagonal(center, dia);
            center->setF({centerForce[0]+f_ij[0], centerForce[1]+f_ij[1], centerForce[2]+f_ij[2]});
            dia->setF({otherForce[0]-f_ij[0], otherForce[1]-f_ij[1], otherForce[2]-f_ij[2]});
        }
        pair_count ++;
    }
};

void Membrane::setMovingParticles(std::vector<Particle*> movPart){
	movingParticles = movPart;
}

int Membrane::getMovingParticleCount(){
    return movingParticles.size();
}