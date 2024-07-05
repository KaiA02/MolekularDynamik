//
// Created by jh on 04.07.2024.
//
#include "MembranePair.h"



MembranePair::MembranePair(Particle* p){
    center = p;
    std::vector<Particle*> dir{};
    direct = dir;
    std::vector<Particle*> dia{};
    diagonal = dia;
}

Particle* MembranePair::getCenter(){
    return center;
}

std::vector<Particle*> MembranePair::getDirect(){
    return direct;
}

std::vector<Particle*> MembranePair::getDiagonal(){
    return diagonal;
}

void MembranePair::appendDirect(Particle* p){
    direct.push_back(p);
};

void MembranePair::appendDiagonal(Particle* p){
    diagonal.push_back(p);
};