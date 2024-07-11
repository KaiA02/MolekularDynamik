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
    for(auto pair: pairs){
        Particle* center = pair.getCenter();
        for(auto dir: pair.getDirect()){
            double distance = calcDistance(center->getX(), dir->getX());
            spdlog::warn("distance {}", distance);
        }
        for(auto dia: pair.getDiagonal()){
            double distance = calcDistance(center->getX(), dia->getX());
            spdlog::warn("distance {}", distance);
        }
    }

}
double Membrane::calcDistance(std::array<double, 3> x1, std::array<double, 3> x2){
    return std::sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1])+(x1[2]-x2[2])*(x1[2]-x2[2]));
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
    //getAveragePairSize();
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
    applyMovement();
    std::vector<double> f_ij;
    std::array<double, 3> f_ij_result;
    Particle* center;
    std::array<double,3> centerF;
    double r0 = 2.2;
    double sqr0 = 2.2 * sqrt(2);
    for(auto pair: pairs){
        center = pair.getCenter();
        f_ij_result = center->getF();
        for(auto dir: pair.getDirect()){
            f_ij = calc.calculateHarmonicForce(center, dir, r0);
            f_ij_result = {f_ij_result[0] + f_ij[0], f_ij_result[1] + f_ij[1], f_ij_result[2] + f_ij[2]};
        }
        for(auto dia: pair.getDiagonal()){
            f_ij = calc.calculateHarmonicForce(center, dia, sqr0);
            f_ij_result = {f_ij_result[0] + f_ij[0], f_ij_result[1] + f_ij[1], f_ij_result[2] + f_ij[2]};
        }
        center->setF(f_ij_result);
    }
    //getAveragePairSize();
};

void Membrane::setMovingParticles(std::vector<Particle*> movPart){
	movingParticles = movPart;
}

int Membrane::getMovingParticleCount(){
    return movingParticles.size();
}