//
// Created by kaiarenja on 16.05.24.
//

#include "Calculations.h"

#include <cmath>

#include "Container/ParticleContainer.h"
#include "Particle.h"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"
#include <omp.h>

Calculations::Calculations(BaseParticleContainer &other) : particles(other) {
  r_cutoff = std::numeric_limits<double>::infinity();
  smoothLJ = false;
}

void Calculations::setR_cutoff(double r) { r_cutoff = r; }
void Calculations::setG_grav(double g) { g_grav = g; }

void Calculations::setSmoothLJ(bool SLJ) { smoothLJ = SLJ; }

void Calculations::setR_L(double R_L) { r_l = R_L; }

void Calculations::setParallelStrategy(int strat){ParallelStrategy = strat;}

void Calculations::calculateX(double delta_t) {
    switch (ParallelStrategy) {
        case 0: // No parallel execution
            for (size_t i = 0; i < particles.size(); ++i) {
                auto &p = particles.getParticles()[i];
                std::array<double, 3> X = p.getX();
                std::array<double, 3> V = p.getV();
                std::array<double, 3> F = p.getF();
                double M = p.getM();

                std::array<double, 3> newPosition = {
                    X[0] + delta_t * V[0] + delta_t * delta_t * (F[0] / (2 * M)),
                    X[1] + delta_t * V[1] + delta_t * delta_t * (F[1] / (2 * M)),
                    X[2] + delta_t * V[2] + delta_t * delta_t * (F[2] / (2 * M))
                };

                p.setX(newPosition);
            }
            break;

        case 1: // Parallel execution using OpenMP for loops
            #pragma omp parallel for
            for (size_t i = 0; i < particles.size(); ++i) {
                auto &p = particles.getParticles()[i];
                std::array<double, 3> X = p.getX();
                std::array<double, 3> V = p.getV();
                std::array<double, 3> F = p.getF();
                double M = p.getM();

                std::array<double, 3> newPosition = {
                    X[0] + delta_t * V[0] + delta_t * delta_t * (F[0] / (2 * M)),
                    X[1] + delta_t * V[1] + delta_t * delta_t * (F[1] / (2 * M)),
                    X[2] + delta_t * V[2] + delta_t * delta_t * (F[2] / (2 * M))
                };

                p.setX(newPosition);
            }
            break;

        case 2: // Parallel execution using OpenMP task parallelism
            #pragma omp parallel
            {
                #pragma omp single
                {
                    for (size_t i = 0; i < particles.size(); ++i) {
                        #pragma omp task firstprivate(i)
                        {
                            auto &p = particles.getParticles()[i];
                            std::array<double, 3> X = p.getX();
                            std::array<double, 3> V = p.getV();
                            std::array<double, 3> F = p.getF();
                            double M = p.getM();

                            std::array<double, 3> newPosition = {
                                X[0] + delta_t * V[0] + delta_t * delta_t * (F[0] / (2 * M)),
                                X[1] + delta_t * V[1] + delta_t * delta_t * (F[1] / (2 * M)),
                                X[2] + delta_t * V[2] + delta_t * delta_t * (F[2] / (2 * M))
                            };

                            p.setX(newPosition);
                        }
                    }
                }
            }
            break;

        default:
            throw std::invalid_argument("Unknown parallel strategy");
    }
}

void Calculations::calculateV(double delta_t) {
    switch (ParallelStrategy) {
        case 0: // No parallel execution
            for (size_t i = 0; i < particles.size(); ++i) {
                auto &p = particles.getParticles()[i];
                std::array<double, 3> V = p.getV();
                std::array<double, 3> F = p.getF();
                std::array<double, 3> oldF = p.getOldF();
                double M = p.getM();

                std::array<double, 3> newVelocity = {
                    V[0] + (delta_t / (2 * M)) * (F[0] + oldF[0]),
                    V[1] + (delta_t / (2 * M)) * (F[1] + oldF[1]),
                    V[2] + (delta_t / (2 * M)) * (F[2] + oldF[2])
                };

                p.setV(newVelocity);
            }
            break;

        case 1: // Parallel execution using OpenMP for loops
            #pragma omp parallel for
            for (size_t i = 0; i < particles.size(); ++i) {
                auto &p = particles.getParticles()[i];
                std::array<double, 3> V = p.getV();
                std::array<double, 3> F = p.getF();
                std::array<double, 3> oldF = p.getOldF();
                double M = p.getM();

                std::array<double, 3> newVelocity = {
                    V[0] + (delta_t / (2 * M)) * (F[0] + oldF[0]),
                    V[1] + (delta_t / (2 * M)) * (F[1] + oldF[1]),
                    V[2] + (delta_t / (2 * M)) * (F[2] + oldF[2])
                };

                p.setV(newVelocity);
            }
            break;

        case 2: // Parallel execution using OpenMP task parallelism
            #pragma omp parallel
            {
                #pragma omp single
                {
                    for (size_t i = 0; i < particles.size(); ++i) {
                        #pragma omp task firstprivate(i)
                        {
                            auto &p = particles.getParticles()[i];
                            std::array<double, 3> V = p.getV();
                            std::array<double, 3> F = p.getF();
                            std::array<double, 3> oldF = p.getOldF();
                            double M = p.getM();

                            std::array<double, 3> newVelocity = {
                                V[0] + (delta_t / (2 * M)) * (F[0] + oldF[0]),
                                V[1] + (delta_t / (2 * M)) * (F[1] + oldF[1]),
                                V[2] + (delta_t / (2 * M)) * (F[2] + oldF[2])
                            };

                            p.setV(newVelocity);
                        }
                    }
                }
            }
            break;

        default:
            throw std::invalid_argument("Unknown parallel strategy");
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
        newForce = {prefactor * (x2[0] - x1[0]), prefactor * (x2[1] - x1[1]),
                    prefactor * (x2[2] - x1[2])};
      }
    }
    p1.setF(newForce);
  }
}


void Calculations::LCcalculateLJF(std::vector<Particle*> &center, std::vector<Particle> &other, std::vector<EpsilonSigma> EAndS) {
    if (center.size() > 1) {
        calculateLJFcenter(center, EAndS);
    }

    if (other.size() > 0) {
        switch (ParallelStrategy) {
            case 0: // No parallel execution
                for (size_t i = 0; i < center.size(); ++i) {
                    Particle* pi = center[i];
                    std::array<double, 3> newForcei = {0, 0, 0};

                    for (size_t j = 0; j < other.size(); ++j) {
                        Particle& pj = other[j];
                        std::array<double, 3> f_ij{};
                        double e = pi->getEpsilon();
                        double s = pi->getSigma();

                        if (pi->getType() == 0 && pj.getType() == 0) { // is Membrane
                            if (!pi->isNeighbour(&pj)) {
                                f_ij = decideForceMethod(pi, &pj, e, s);
                            }
                        } else {
                            if (pi->getType() != pj.getType()) {
                                for (auto& entry : EAndS) {
                                    if (entry.isRight(pi->getType(), pj.getType())) {
                                        e = entry.getEpsilon();
                                        s = entry.getSigma();
                                        break;
                                    }
                                }
                            }
                            f_ij = decideForceMethod(pi, &pj, e, s);
                        }

                        newForcei[0] += f_ij[0];
                        newForcei[1] += f_ij[1];
                        newForcei[2] += f_ij[2];
                    }

                    auto force = pi->getF();
                    force[0] += newForcei[0];
                    force[1] += newForcei[1];
                    force[2] += newForcei[2];
                    pi->setF(force);
                }
                break;

            case 1: // Parallel execution using OpenMP for loops
                #pragma omp parallel
                {
                    #pragma omp for nowait
                    for (size_t i = 0; i < center.size(); ++i) {
                        Particle* pi = center[i];
                        std::array<double, 3> newForcei = {0, 0, 0};

                        for (size_t j = 0; j < other.size(); ++j) {
                            Particle& pj = other[j];
                            std::array<double, 3> f_ij{};
                            double e = pi->getEpsilon();
                            double s = pi->getSigma();

                            if (pi->getType() == 0 && pj.getType() == 0) { // is Membrane
                                if (!pi->isNeighbour(&pj)) {
                                    f_ij = decideForceMethod(pi, &pj, e, s);
                                }
                            } else {
                                if (pi->getType() != pj.getType()) {
                                    for (auto& entry : EAndS) {
                                        if (entry.isRight(pi->getType(), pj.getType())) {
                                            e = entry.getEpsilon();
                                            s = entry.getSigma();
                                            break;
                                        }
                                    }
                                }
                                f_ij = decideForceMethod(pi, &pj, e, s);
                            }

                            newForcei[0] += f_ij[0];
                            newForcei[1] += f_ij[1];
                            newForcei[2] += f_ij[2];
                        }

                        #pragma omp critical
                        {
                            auto force = pi->getF();
                            force[0] += newForcei[0];
                            force[1] += newForcei[1];
                            force[2] += newForcei[2];
                            pi->setF(force);
                        }
                    }
                }
                break;

            case 2: // Parallel execution using OpenMP task parallelism
                #pragma omp parallel
                {
                    #pragma omp single
                    {
                        for (size_t i = 0; i < center.size(); ++i) {
                            #pragma omp task firstprivate(i)
                            {
                                Particle* pi = center[i];
                                std::array<double, 3> newForcei = {0, 0, 0};

                                for (size_t j = 0; j < other.size(); ++j) {
                                    Particle& pj = other[j];
                                    std::array<double, 3> f_ij{};
                                    double e = pi->getEpsilon();
                                    double s = pi->getSigma();

                                    if (pi->getType() == 0 && pj.getType() == 0) { // is Membrane
                                        if (!pi->isNeighbour(&pj)) {
                                            f_ij = decideForceMethod(pi, &pj, e, s);
                                        }
                                    } else {
                                        if (pi->getType() != pj.getType()) {
                                            for (auto& entry : EAndS) {
                                                if (entry.isRight(pi->getType(), pj.getType())) {
                                                    e = entry.getEpsilon();
                                                    s = entry.getSigma();
                                                    break;
                                                }
                                            }
                                        }
                                        f_ij = decideForceMethod(pi, &pj, e, s);
                                    }

                                    newForcei[0] += f_ij[0];
                                    newForcei[1] += f_ij[1];
                                    newForcei[2] += f_ij[2];
                                }

                                #pragma omp critical
                                {
                                    auto force = pi->getF();
                                    force[0] += newForcei[0];
                                    force[1] += newForcei[1];
                                    force[2] += newForcei[2];
                                    pi->setF(force);
                                }
                            }
                        }
                    }
                }
                break;

            default:
                throw std::invalid_argument("Unknown parallel strategy");
        }
    }
}






void Calculations::calculateLJFcenter(std::vector<Particle *> &center,
                                      const std::vector<EpsilonSigma> EAndS) {

  std::array<double, 3> f_ij{};
  double e = 1;
  double s = 1;
  Particle* pi;
  Particle* pk;
  for (size_t i = 0; i < center.size() - 1; ++i) {
    pi = center.at(i);
    for (size_t j = i + 1; j < center.size(); ++j) {
      pk = center.at(j);
      if(pi->getType() == 0 && pk->getType() == 0) {
        if(!pi->isNeighbour(pk)){
          f_ij = decideForceMethod(pi, pk, pi->getEpsilon(), pi->getSigma());
          pi->setF(pi->getF() + f_ij);
          pk->setF(pk->getF() - f_ij);
          spdlog::debug("calculated Membrane LJF {} {} {}", f_ij[0], f_ij[1], f_ij[2]);
		}
      } else {
        if(pi->getType() != pk->getType()) {
          for(auto entry : EAndS) {
            if(entry.isRight(pi->getType(), pk->getType())) {
              e = entry.getEpsilon();
              s = entry.getSigma();
              break;
            }
          }
        } else {
          e = pi->getEpsilon();
          s = pi->getSigma();
        }
        f_ij = decideForceMethod(pi, pk, e, s);
        pi->setF(pi->getF() + f_ij);
        pk->setF(pk->getF() - f_ij);
      }
    }
  }
}



std::array<double, 3> Calculations::calculateLJF(Particle *p1, Particle *p2, double e, double s) {
  std::array<double, 3> x1 = p1->getX();
  std::array<double, 3> x2 = p2->getX();
  std::array<double, 3> f_ij{};
  std::array<double, 3> displacement_vector = { x1[0] - x2[0], x1[1] - x2[1], x1[2] - x2[2] };

  double distance_squared = displacement_vector[0] * displacement_vector[0] +
                            displacement_vector[1] * displacement_vector[1] +
                            displacement_vector[2] * displacement_vector[2];
  double distance = std::sqrt(distance_squared);

  if (distance <= r_cutoff) {
    double s_over_d = s / distance;
    double s_over_d2 = s_over_d * s_over_d;
    double s_over_d6 = s_over_d2 * s_over_d2 * s_over_d2;
    double s_over_d12 = s_over_d6 * s_over_d6;

    double forcefactor = (-24 * e / distance_squared) * (s_over_d6 - 2 * s_over_d12);

    f_ij = { forcefactor * displacement_vector[0],
             forcefactor * displacement_vector[1],
             forcefactor * displacement_vector[2] };
  }

  return f_ij;
}

std::array<double, 3> Calculations::calculateSmoothLJF(Particle *p1,
                                                       Particle *p2, double e,
                                                       double s) {
  std::array<double, 3> x1 = p1->getX();
  std::array<double, 3> x2 = p2->getX();
  std::array<double, 3> f_ij{};
  std::array<double, 3> displacement_vector = {x1[0] - x2[0], x1[1] - x2[1],
                                               x1[2] - x2[2]};
  double distance = sqrt(displacement_vector[0] * displacement_vector[0] +
                         displacement_vector[1] * displacement_vector[1] +
                         displacement_vector[2] * displacement_vector[2]);

  double distance_exp2 = distance * distance;
  double distance_exp6 = distance_exp2 * distance_exp2 * distance_exp2;
  double distance_exp14 = distance_exp6 * distance_exp6 * distance_exp2;
  double sigma_exp6 = pow(s, 6);

  double forcefactor1 =
      -(24 * sigma_exp6 * e / (distance_exp14 * pow(r_cutoff - r_l, 3)) *
        (r_cutoff - distance));

  double forcefactor2_1 =
      r_cutoff * r_cutoff * (2 * sigma_exp6 - distance_exp6) +
      r_cutoff * (3 * r_l - distance) * (distance_exp6 - 2 * sigma_exp6);

  double forcefactor2_2_1 = 5 * r_l * sigma_exp6 - 2 * r_l * distance_exp6 -
                            3 * sigma_exp6 * distance +
                            distance_exp6 * distance;

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
        double distance = calcDistance(x1, x2);
        if(distance <= r_l) {
            return calculateLJF(p1, p2, e, s);
        } else if(r_l <= distance && distance <= r_cutoff) {
            return calculateSmoothLJF(p1, p2, e, s);
        } else {
            return {0.0,0.0,0.0};
        }
    }else {
        return calculateLJF(p1, p2, e, s);
    }
}

std::vector<double> Calculations::calculateHarmonicForce(Particle *p1, Particle *p2, double r0) {
  std::array<double, 3> xi = p1->getX();
  std::array<double, 3> xj = p2->getX();
  std::vector<double> i_j = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
  double distance = calcDistance(xi, xj);
  double scalar = stiffness * (distance - r0);
  for (int i = 0; i < 3; ++i) {
    i_j[i] = scalar * (xj[i] - xi[i]) / distance;
  }
  return i_j;
}



double Calculations::calcDistance(std::array<double, 3> x1, std::array<double, 3> x2){
  double dx0 = x1[0] - x2[0];
  double dx1 = x1[1] - x2[1];
  double dx2 = x1[2] - x2[2];
  return std::sqrt(dx0 * dx0 + dx1 * dx1 + dx2 * dx2);

}

double Calculations::calculateDistanceBetweenParticles(Particle *p1, Particle *p2) {
  std::array<double, 3> x1 = p1->getX();
  std::array<double, 3> x2 = p2->getX();
  double dx0 = x1[0] - x2[0];
  double dx1 = x1[1] - x2[1];
  double dx2 = x1[2] - x2[2];
  return std::sqrt(dx0 * dx0 + dx1 * dx1 + dx2 * dx2);
}

double Calculations::calculateDiffusion(std::vector<Particle> particles,std::vector<Particle> prevParticles) {
  double diffusion = 0;
  for (int i = 0; i < particles.size(); i++) {
    Particle currentParticle = particles.at(i);
    Particle prevParticle = prevParticles.at(i);
    diffusion +=
        calculateDistanceBetweenParticles(&currentParticle, &prevParticle);
    ;
  }
  return diffusion / particles.size();
}

std::vector<std::vector<double>> Calculations::computeDistances(std::vector<Particle> particles) {
  int numParticles = particles.size();

  // Initialize a 2D vector to store distances
  std::vector<std::vector<double>> distances(
      numParticles, std::vector<double>(numParticles, 0.0));

  // Calculate pairwise distances and store them
  for (int i = 0; i < numParticles; ++i) {
    for (int j = i + 1; j < numParticles; ++j) {
      double distance =
          calculateDistanceBetweenParticles(&particles[i], &particles[j]);
      distances[i][j] = distance;
      distances[j][i] = distance; // Symmetric matrix
      if (distance > maxDistance) {
        maxDistance = distance;
      }
    }
  }
  return distances;
}

std::map<double, double> Calculations::calculateLocalDensities(const std::vector<Particle> particles, double deltaR) {
  std::vector<std::vector<double>> distances = computeDistances(particles);
  std::map<double, double> localDensities;
  int numIntervals = static_cast<int>(std::ceil(maxDistance / deltaR));
  std::vector<int> intervalCounts(numIntervals, 0);
  double currentDistance = 0.0;

  for (size_t i = 0; i < distances.size(); ++i) {
    for (size_t j = i + 1; j < distances[i].size(); ++j) {
      currentDistance = distances[i][j];
      int intervalIndex = static_cast<int>(currentDistance / deltaR);
      if (intervalIndex < numIntervals) {
        intervalCounts[intervalIndex]++;
      }
    }
  }

  double ri = 0.0;
  for (int m = 0; m < numIntervals; ++m) {
    ri = m * deltaR;
    double volume = (4.0 / 3.0) * M_PI * (pow(ri + deltaR, 3) - pow(ri, 3));
    double rdf = static_cast<double>(intervalCounts[m]) / volume;

    // store the local density for interval ri
    localDensities[ri] = rdf;
  }

  return localDensities;

}