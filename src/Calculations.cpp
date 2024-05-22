//
// Created by kaiarenja on 16.05.24.
//

#include "Calculations.h"

#include <cmath>

#include "Particle.h"
#include "ParticleContainer.h"

Calculations::Calculations(ParticleContainer &other) : particles(other) {}

ParticleContainer &Calculations::getParticles() { return particles; }

/**
 * @brief the following function calculates the new positions
 *
 * therefore we use a for-loops
 * to calculate the physics behind the positions (Velocity-Störmer-Verlet)
 */
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

/**
 * @brief the following function calculates the new velocities
 *
 * therefore we use a for-loops
 * to calculate the physics behind the velocities (Velocity-Störmer-Verlet)
 */
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

/**
 * @brief the following function calculates the gravitational force
 *
 * therefore we use an interator and also for-loops
 * to calculate the physics behind the forces
 */
void Calculations::calculateF() {
  for (auto &p1 : particles) {
    std::array<double, 3> newForce = {0.0, 0.0, 0.0};

    for (auto &p2 : particles) {
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

/**
 * @brief the following function calculates the Leonard Jones force
 */
void Calculations::calculateLJF() {
  for (int i = 0; i < particles.size(); ++i) {
    Particle &pi = particles.getParticles().at(i);
    pi.setF({0,0,0});
  }
  for (int i = 0; i < particles.size() - 1; ++i) {
    Particle &pi = particles.getParticles().at(i);
    for (int j = i + 1; j < particles.size(); ++j) {
      Particle &pj = particles.getParticles().at(j);
      std::array<double, 3> displacement_vector = {
          pi.getX().at(0) - pj.getX().at(0), pi.getX().at(1) - pj.getX().at(1),
          pi.getX().at(2) - pj.getX().at(2)};
      double distance = sqrt(pow(displacement_vector.at(0), 2) +
                             pow(displacement_vector.at(1), 2) +
                             pow(displacement_vector.at(2), 2));

      std::array<double, 3> f_ij = {
          -1 * ((24 * epsilion) / (pow(distance, 2))) *
              (pow((sigma) / distance, 6) - 2 * pow((sigma) / distance, 12)) *
              displacement_vector.at(0),
          -1 * ((24 * epsilion) / (pow(distance, 2))) *
              (pow((sigma) / distance, 6) - 2 * pow((sigma) / distance, 12)) *
              displacement_vector.at(1),
          -1 * ((24 * epsilion) / (pow(distance, 2))) *
              (pow((sigma) / distance, 6) - 2 * pow((sigma) / distance, 12)) *
              displacement_vector.at(2)};
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
