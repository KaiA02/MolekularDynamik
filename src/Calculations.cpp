//
// Created by kaiarenja on 16.05.24.
//

#include "Calculations.h"

#include <cmath>

#include "Particle.h"
#include "Container/ParticleContainer.h"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

Calculations::Calculations(BaseParticleContainer &other) : particles(other) {
  r_cutoff = std::numeric_limits<double>::infinity();
}

void Calculations::setR_cutoff(double r) {
  r_cutoff = r;
}
 void Calculations::setG_grav(double g) {
   g_grav = g;
 }

void Calculations::setSmoothLJ(bool SLJ) {
  smoothLJ = SLJ;
}

void Calculations::setR_L(double R_L) {
    r_l = R_L;
}


void Calculations::calculateX(double delta_t) {
  std::array<double, 3>  newPosition{};
  std::array<double, 3> X{};
  std::array<double, 3> V{};
  std::array<double, 3> F{};
  double M = 0;

  for (auto &p : particles) {
    X = p.getX();
    V = p.getV();
    F = p.getF();
    M = p.getM();

    newPosition = { X[0] + delta_t * V[0] + delta_t * delta_t * (F[0] / (2 * M)),
                      X[1] + delta_t * V[1] + delta_t * delta_t * (F[1] / (2 * M)),
                      X[2] + delta_t * V[2] + delta_t * delta_t * (F[2] / (2 * M))};
    p.setX(newPosition);
  }
}


void Calculations::calculateV(double delta_t) {
  std::array<double, 3> newVelocity{};
  std::array<double, 3> V{};
  std::array<double, 3> F{};
  std::array<double, 3> oldF{};
  double M = 0;
  for (auto &p : particles) {
    V = p.getV();
    F = p.getF();
    oldF = p.getOldF();
    M = p.getM();

    newVelocity = { V[0] + delta_t * (F[0] + oldF[0]) / (2 * M),
                      V[1] + delta_t * (F[1] + oldF[1]) / (2 * M),
                      V[2] + delta_t * (F[2] + oldF[2]) / (2 * M)};
    p.setV(newVelocity);
  }
}

void Calculations::calculateF() {
  std::array<double, 3> newForce;
  double distCubed;
  double prefactor;
  std::array<double, 3> x1;
  std::array<double, 3> x2;
  for (int k = 0; k < particles.size(); k++) {
    Particle &p1 = particles.getParticles().at(k);
    newForce = {0.0, 0.0, 0.0};
    x1 = p1.getX();
    for (int j = 0; j < particles.size(); j++) {
      Particle &p2 = particles.getParticles().at(j);
      x2 = p2.getX();
      if (&p1 != &p2) {
        distCubed = std::pow(ArrayUtils::L2Norm(x1 - x2), 3);
        prefactor = (p1.getM() * p2.getM()) / distCubed;
        newForce = {prefactor * (x2[0] - x1[0]),
                      prefactor * (x2[1] - x1[1]),
                      prefactor * (x2[2] - x1[2])};
      }
    }
    p1.setF(newForce);
  }
}



void Calculations::LCcalculateLJF(std::vector<Particle*> &center, std::vector<Particle> &other, std::vector<EpsilonSigma> EAndS) {

  std::array<double, 3> f_ij;
  std::array<double, 3> newForcei;
  double e;
  double s;

  if (center.size() > 1) {
    calculateLJFcenter(center, EAndS);
  }
  if (other.size() > 0) {
    for(auto pi : center){
      newForcei = pi->getF();
      for(auto pj : other) {
        if(pi->getType() != pj.getType()) {
          for(auto entry : EAndS) {
          if(entry.isRight(pi->getType(), pj.getType())) {
            e = entry.getEpsilon();
            s = entry.getSigma();
            break;
          }
        }
        } else {
          e = pi->getEpsilon();
          s = pi->getSigma();
        }
        f_ij = decideForceMethod(pi, &pj, e, s);

        newForcei = {newForcei.at(0) + f_ij.at(0), newForcei.at(1) + f_ij.at(1), newForcei.at(2) + f_ij.at(2)};
      }
      pi->setF(newForcei);
    }
  }
}


void Calculations::calculateLJFcenter(std::vector<Particle *> &center, const std::vector<EpsilonSigma> EAndS) {

  std::array<double, 3> f_ij{};
  double e = 1;
  double s = 1;

  for (size_t i = 0; i < center.size() - 1; ++i) {
    Particle *pi = center.at(i);
    for (size_t j = i + 1; j < center.size(); ++j) {
      Particle *pj = center.at(j);
      if(pi->getType() != pj->getType()) {
        for(auto entry : EAndS) {
          if(entry.isRight(pi->getType(), pj->getType())) {
            e = entry.getEpsilon();
            s = entry.getSigma();
            break;
          }
        }
      } else {
        e = pi->getEpsilon();
        s = pi->getSigma();
      }
      f_ij = decideForceMethod(pi, pj, e, s);
      pi->setF(pi->getF() + f_ij);
      pj->setF(pj->getF() - f_ij);
    }
  }
}

std::array<double, 3> Calculations::calculateLJF(Particle *p1, Particle *p2, double e, double s) {
  std::array<double,3> x1 = p1->getX();
  std::array<double,3> x2 = p2->getX();
  std::array<double, 3> f_ij{};
  std::array<double, 3> displacement_vector = { x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]};
  double distance = sqrt(displacement_vector[0] * displacement_vector[0]+
                         displacement_vector[1] * displacement_vector[1] +
                         displacement_vector[2] * displacement_vector[2]);
if (distance <= r_cutoff) {
  double forcefactor =
      ((-24 * e) / pow(distance, 2)) *
      (pow((s) / distance, 6) - 2 * pow((s) / distance, 12));

  f_ij = {forcefactor * displacement_vector[0],
          forcefactor * displacement_vector[1],
          forcefactor * displacement_vector[2]};

}
  return f_ij;
}


std::array<double, 3> Calculations::calculateSmoothLJF(Particle *p1, Particle *p2, double e, double s) {
    std::array<double,3> x1 = p1->getX();
    std::array<double,3> x2 = p2->getX();
    std::array<double, 3> f_ij{};
    std::array<double, 3> displacement_vector = { x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]};
    double distance = sqrt(displacement_vector[0] * displacement_vector[0]+
                           displacement_vector[1] * displacement_vector[1] +
                           displacement_vector[2] * displacement_vector[2]);

    double distance_exp2 = distance * distance;
    double distance_exp6 = distance_exp2 * distance_exp2 * distance_exp2;
    double distance_exp14 = distance_exp6 * distance_exp6 * distance_exp2;
    double sigma_exp6 = pow(s, 6);

    double forcefactor1 = -(24 * sigma_exp6 * e / (distance_exp14 * pow(r_cutoff - r_l, 3)) * (r_cutoff - distance));

    double forcefactor2_1 = r_cutoff * r_cutoff * (2 * sigma_exp6 - distance_exp6) +
                            r_cutoff * (3 * r_l - distance) * (distance_exp6 - 2 * sigma_exp6);

    double forcefactor2_2_1 = 5 * r_l * sigma_exp6 - 2 * r_l * distance_exp6 -
                              3 * sigma_exp6 * distance + distance_exp6 * distance;

    double forcefactor2_2 = distance * forcefactor2_2_1;

    double forcefactor = forcefactor1 * (forcefactor2_1 + forcefactor2_2);

    f_ij = {forcefactor * -displacement_vector[0],
            forcefactor * -displacement_vector[1],
            forcefactor * -displacement_vector[2]};

    return f_ij;
}

std::array<double, 3> Calculations::decideForceMethod(Particle *p1, Particle *p2, double e, double s) {
    if(smoothLJ) {
        std::array<double,3> x1 = p1->getX();
        std::array<double,3> x2 = p2->getX();
        std::array<double, 3> displacement_vector = { x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]};
        double distance = sqrt(displacement_vector[0] * displacement_vector[0]+
                               displacement_vector[1] * displacement_vector[1] +
                               displacement_vector[2] * displacement_vector[2]);
        if(distance <= r_l) {
            return calculateLJF(p1, p2, e, s);
        } else if(r_l <= distance && distance <= r_cutoff) {
            return calculateSmoothLJF(p1, p2, e, s);
        } else {
            return {};
        }
    }else {
        return calculateLJF(p1, p2, e, s);
    }
}

std::array<double, 3> Calculations::calculateHarmonicForce(Particle *p1, Particle *p2){
  std::array<double,3> x1 = p1->getX();
  std::array<double,3> x2 = p2->getX();
  std::array<double, 3> f_ij{};
  std::array<double, 3> displacement_vector = { x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]};
  double distance = sqrt(displacement_vector[0] * displacement_vector[0]+
                         displacement_vector[1] * displacement_vector[1] +
                         displacement_vector[2] * displacement_vector[2]);
  double forcefactor = stiffness * (distance - avgBondLength);
  f_ij = {forcefactor * (-displacement_vector[0] / distance),
          forcefactor * (-displacement_vector[1] / distance),
          forcefactor * (-displacement_vector[2] / distance)};

  return f_ij;

}

std::array<double, 3> Calculations::calculateHarmonicForceDiagonal(Particle *p1, Particle *p2){
  std::array<double,3> x1 = p1->getX();
  std::array<double,3> x2 = p2->getX();
  std::array<double, 3> f_ij{};
  std::array<double, 3> displacement_vector = { x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2]};
  double distance = sqrt(displacement_vector[0] * displacement_vector[0]+
                         displacement_vector[1] * displacement_vector[1] +
                         displacement_vector[2] * displacement_vector[2]);
  double forcefactor = stiffness * (distance - (sqrt(2) * avgBondLength));
  f_ij = {forcefactor * (-displacement_vector[0] / distance),
          forcefactor * (-displacement_vector[1] / distance),
          forcefactor * (-displacement_vector[2] / distance)};
  return f_ij;
}