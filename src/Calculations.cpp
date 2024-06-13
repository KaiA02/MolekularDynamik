//
// Created by kaiarenja on 16.05.24.
//

#include "Calculations.h"

#include <cmath>

#include "Particle.h"
#include "Container/ParticleContainer.h"
#include "spdlog/spdlog.h"


Calculations::Calculations(ParticleContainer &other) : particles(other) {}

Calculations::Calculations(LCParticleContainer *other): particles(*other) {
}


void Calculations::calculateX(double delta_t) {
  for (auto &p : particles) {
    std::array<double, 3> newPosition;
    for (int i = 0; i < 3; ++i) {
      newPosition[i] = p.getX()[i] + delta_t * p.getV()[i] +
                       delta_t * delta_t * (p.getF()[i] / (2 * p.getM()));
    }
    p.setX(newPosition);
  }
}


void Calculations::calculateV(double delta_t) {
  for (auto &p : particles) {
    std::array<double, 3> newVelocity;
    for (int i = 0; i < 3; ++i) {
      newVelocity[i] = p.getV()[i] + delta_t * (p.getF()[i] + p.getOldF()[i]) /
                                         (2 * p.getM());
    }
    p.setV(newVelocity);
  }
}

void Calculations::calculateF() {
  for (int k = 0; k < particles.size(); k++) {
    Particle &p1 = particles.getParticles().at(k);
    std::array<double, 3> newForce = {0.0, 0.0, 0.0};
    for (int j = 0; j < particles.size(); j++) {
      Particle &p2 = particles.getParticles().at(j);
      if (&p1 != &p2) {
        double distSquared = 0.0;
        for (int i = 0; i < 3; ++i) {
          distSquared += std::pow(p1.getX()[i] - p2.getX()[i], 2);
        }

        double distCubed = std::pow(distSquared, 1.5);

        double prefactor = (p1.getM() * p2.getM()) / distCubed;

        for (int i = 0; i < 3; ++i) {
          newForce[i] += prefactor * (p2.getX()[i] - p1.getX()[i]);
        }
      }
    }
    p1.setF(newForce);

  }
}

void Calculations::calculateLJF() {
  for (int i = 0; i < particles.size(); ++i) {
    Particle &pi = particles.getParticles().at(i);
    pi.setF({0, 0, 0});
  }
  for (int i = 0; i < particles.size() - 1; ++i) {
    Particle &pi = particles.getParticles().at(i);
    for (int j = i + 1; j < particles.size(); ++j) {
      Particle &pj = particles.getParticles().at(j);
      double fake_r_cutoff = 10000000000000000;
      std::array<double, 3> f_ij = calculateLJF(&pi, &pj, fake_r_cutoff);

      std::array<double, 3> newForcei = {0.0, 0.0, 0.0};
      std::array<double, 3> newForcej = {0.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k) {
        newForcei.at(k) = pi.getF().at(k) + f_ij.at(k);
        newForcej.at(k) = pj.getF().at(k) - f_ij.at(k);
      }
      pi.setF(newForcei);
      pj.setF(newForcej);
    }
  }
}


void Calculations::LCcalculateLJF(std::vector<Particle*> &center, std::vector<Particle*> &other, double r_cutoff) {
  if (center.size() > 1) {
    calculateLJFcenter(center, r_cutoff);
  }
  if (other.size() > 0) {
   for(auto pi : center) {
      for(auto pj : center) {
        std::array<double, 3> f_ij = calculateLJF(pi, pj, r_cutoff);
        std::array<double, 3> newForcei = {0.0, 0.0, 0.0};
        //std::array<double, 3> newForcej = {0.0, 0.0, 0.0};
        for (int k = 0; k < 3; ++k) {
          newForcei.at(k) = pi->getF().at(k) + f_ij.at(k);
          // newForcej.at(k) = pj.getF().at(k) - (f_ij.at(k)/2);
        }
        pi->setF(newForcei);
        // pj.setF(newForcej);
      }
    }
  }
}


void Calculations::calculateLJFcenter(std::vector<Particle *> &center, double r_cutoff) {
  for (size_t i = 0; i < center.size() - 1; ++i) {
    Particle *pi = center.at(i);
    for (size_t j = i + 1; j < center.size(); ++j) {
      Particle *pj = center.at(j);
      std::array<double, 3> f_ij = calculateLJF(pi, pj, r_cutoff);
      std::array<double, 3> newForcei = {0.0, 0.0, 0.0};
      std::array<double, 3> newForcej = {0.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k) {
        newForcei.at(k) = pi->getF().at(k) + f_ij.at(k);
        newForcej.at(k) = pj->getF().at(k) - f_ij.at(k);
      }
      pi->setF(newForcei);
      pj->setF(newForcej);
    }
  }
}

std::array<double, 3> Calculations::calculateLJF(Particle *p1, Particle *p2, double r_cutoff) {
  std::array<double, 3> displacement_vector = {
      p1->getX().at(0) - p2->getX().at(0), p1->getX().at(1) - p2->getX().at(1),
      p1->getX().at(2) - p2->getX().at(2)};
  double distance = sqrt(pow(displacement_vector.at(0), 2) +
                         pow(displacement_vector.at(1), 2) +
                         pow(displacement_vector.at(2), 2));
  if(distance <= r_cutoff) {
   double forcefactor =
      ((-24 * epsilion) / pow(distance, 2)) *
      (pow((sigma) / distance, 6) - 2 * pow((sigma) / distance, 12));

  std::array<double, 3> f_ij = {forcefactor * displacement_vector.at(0),
                                forcefactor * displacement_vector.at(1),
                                forcefactor * displacement_vector.at(2)};
  return f_ij;
  } else {
    return {0.0,0.0,0.0};
  }
}
