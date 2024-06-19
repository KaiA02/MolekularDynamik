//
// Created by kaiarenja on 16.05.24.
//

#include "Calculations.h"

#include <cmath>

#include "Particle.h"
#include "Container/ParticleContainer.h"
#include "spdlog/spdlog.h"

Calculations::Calculations(BaseParticleContainer &other) : particles(other) {
  r_cutoff = std::numeric_limits<double>::infinity();
}

void Calculations::setR_cutoff(double r) {
  r_cutoff = r;
}
 void Calculations::setG_grav(double g) {
   g_grav = g;
 }



void Calculations::calculateX(double delta_t) {
  std::array<double, 3> newPosition;
  for (auto &p : particles) {
    newPosition = { p.getX()[0] + delta_t * p.getV()[0] + delta_t * delta_t * (p.getF()[0] / (2 * p.getM())),
                      p.getX()[1] + delta_t * p.getV()[1] + delta_t * delta_t * (p.getF()[1] / (2 * p.getM())),
                      p.getX()[2] + delta_t * p.getV()[2] + delta_t * delta_t * (p.getF()[2] / (2 * p.getM()))};
    p.setX(newPosition);
  }
}


void Calculations::calculateV(double delta_t) {
  std::array<double, 3> newVelocity;
  for (auto &p : particles) {
    newVelocity = { p.getV()[0] + delta_t * (p.getF()[0] + p.getOldF()[0]) / (2 * p.getM()),
                      p.getV()[1] + delta_t * (p.getF()[1] + p.getOldF()[1]) / (2 * p.getM()),
                      p.getV()[2] + delta_t * (p.getF()[2] + p.getOldF()[2]) / (2 * p.getM())};
    p.setV(newVelocity);
  }
}

void Calculations::calculateF() {
  std::array<double, 3> newForce;
  double distSquared;
  double distCubed;
  double prefactor;
  for (int k = 0; k < particles.size(); k++) {
    Particle &p1 = particles.getParticles().at(k);
    newForce = {0.0, 0.0, 0.0};
    for (int j = 0; j < particles.size(); j++) {
      Particle &p2 = particles.getParticles().at(j);
      if (&p1 != &p2) {
        distSquared = std::pow(p1.getX()[0] - p2.getX()[0], 2) + std::pow(p1.getX()[1] - p2.getX()[1], 2) + std::pow(p1.getX()[2] - p2.getX()[2], 2);
        distCubed = std::pow(distSquared, 1.5);
        prefactor = (p1.getM() * p2.getM()) / distCubed;
        newForce = {prefactor * (p2.getX()[0] - p1.getX()[0]),
                      prefactor * (p2.getX()[1] - p1.getX()[1]),
                      prefactor * (p2.getX()[2] - p1.getX()[2])};
      }
    }
    p1.setF(newForce);
  }
}

void Calculations::calculateLJF() {

  std::array<double, 3> f_ij;
  std::array<double, 3> newForcei;
  std::array<double, 3> newForcej;

  for (int i = 0; i < particles.size(); ++i) {
    Particle &pi = particles.getParticles().at(i);
    pi.setF({0, 0, 0});
  }

  for (int i = 0; i < particles.size() - 1; ++i) {
    Particle &pi = particles.getParticles().at(i);
    for (int j = i + 1; j < particles.size(); ++j) {
      Particle &pj = particles.getParticles().at(j);
      f_ij = calculateLJF(&pi, &pj);

      newForcei = {0.0, 0.0, 0.0};
      newForcej = {0.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k) {
        newForcei.at(k) = pi.getF().at(k) + f_ij.at(k);
        newForcej.at(k) = pj.getF().at(k) - f_ij.at(k);
      }
      pi.setF(newForcei);
      pj.setF(newForcej);
    }
  }
}


void Calculations::LCcalculateLJF(std::vector<Particle*> &center, std::vector<Particle> &other) {

  std::array<double, 3> f_ij;
  std::array<double, 3> newForcei;

  if (center.size() > 1) {
    calculateLJFcenter(center);
  }
  if (other.size() > 0) {
    for(auto pi : center){
      newForcei = pi->getF();
      for(auto pj : other) {
        f_ij = calculateLJF(pi, &pj);
        newForcei = {newForcei.at(0) + f_ij.at(0), newForcei.at(1) + f_ij.at(1), newForcei.at(2) + f_ij.at(2)};
      }
      pi->setF(newForcei);
    }
  }
}


void Calculations::calculateLJFcenter(std::vector<Particle *> &center) {

  std::array<double, 3> f_ij;
  std::array<double, 3> newForcei;
  std::array<double, 3> newForcej;

  for (size_t i = 0; i < center.size() - 1; ++i) {
    Particle *pi = center.at(i);
    for (size_t j = i + 1; j < center.size(); ++j) {
      Particle *pj = center.at(j);
      f_ij = calculateLJF(pi, pj);
      newForcei = pi->getF();
      newForcej = pj->getF();
      newForcei = {newForcei.at(0) + f_ij.at(0), newForcei.at(1) + f_ij.at(1), newForcei.at(2) + f_ij.at(2)};
      newForcej = {newForcej.at(0) - f_ij.at(0), newForcej.at(1) - f_ij.at(1), newForcej.at(2) - f_ij.at(2)};
      pi->setF(newForcei);
      pj->setF(newForcej);
    }
  }
}

std::array<double, 3> Calculations::calculateLJF(Particle *p1, Particle *p2) {
  double sigma = (p1->getSigma() + p2->getSigma()) * 0.5;
  double epsilon = std::sqrt(p1->getEpsilon() * p2->getEpsilon());

  std::array<double, 3> displacement_vector = { p1->getX().at(0) - p2->getX().at(0),
                                                p1->getX().at(1) - p2->getX().at(1),
                                                p1->getX().at(2) - p2->getX().at(2)};
  double distance = sqrt(displacement_vector.at(0) * displacement_vector.at(0) +
                         displacement_vector.at(1) * displacement_vector.at(1) +
                         displacement_vector.at(2) * displacement_vector.at(2));
if (distance <= r_cutoff) {
  double forcefactor =
      ((-24 * epsilon) / pow(distance, 2)) *
      (pow((sigma) / distance, 6) - 2 * pow((sigma) / distance, 12));

  std::array<double, 3> f_ij = {forcefactor * displacement_vector.at(0),
                                forcefactor * displacement_vector.at(1),
                                forcefactor * displacement_vector.at(2)};

  return f_ij;

} else {
  return {0.0,0.0,0.0};

}
}
