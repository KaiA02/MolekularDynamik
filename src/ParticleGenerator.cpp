//
// Created by jh on 13.05.2024.
//

#include "ParticleGenerator.h"

#include <stdexcept>

#include "utils//MaxwellBoltzmannDistribution.h"

#include <vector>

ParticleGenerator::ParticleGenerator() {}

void ParticleGenerator::generateCuboid(const Particle &start, int n1, int n2,
                                       int n3, double distance,
                                       double meanVelocity, int dimension, double temp_init) {

    //maxwellVelocity and factor of temperature
    double factor = std::sqrt(temp_init / start.getM());
    std::array<double, 3> addedVelocity = maxwellBoltzmannDistributedVelocity(meanVelocity, dimension);
    addedVelocity = {addedVelocity.at(2) * factor + start.getV()[0], addedVelocity.at(0) * factor + start.getV()[1], addedVelocity.at(1) * factor + start.getV()[2]};

    for (int x = 0; x < n1; x++) {
        for (int y = 0; y < n2; y++) {
            for (int z = 0; z < n3; z++) {
                Particle p({start.getX()[0] + x * distance,
                    start.getX()[1] + y * distance,
                    start.getX()[2] + z * distance},
                   addedVelocity, start.getM(), start.getType(),
                   start.getEpsilon(), start.getSigma());
                allParticles.push_back(p);
      }
    }
  }
}

std::vector<Particle>& ParticleGenerator::getAllParticles() { return allParticles; }

std::vector<Particle*> ParticleGenerator::getAllParticlesPointers(){
std::vector<Particle*> result{};
for(auto& p: allParticles){
    result.push_back(&p);
}
return result;
}


void ParticleGenerator::generateDisk(const Particle &center, int radius, double distance,
                                       double meanVelocity, int dimension, double temp_init) {

    //maxwellVelocity and factor of temperature
    double factor = std::sqrt(temp_init / center.getM());
    std::array<double, 3> addedVelocity = maxwellBoltzmannDistributedVelocity(meanVelocity, dimension);
    addedVelocity = {addedVelocity.at(1) * factor + center.getV()[0], addedVelocity.at(0) * factor + center.getV()[1], addedVelocity.at(2) * factor + center.getV()[2]};

    //maxWellBoltzmann returns huge value for x for no reason. We switch this value to y to support Gravitat
    double x_velocity = addedVelocity.at(0);
    addedVelocity.at(0) = addedVelocity.at(1);
    addedVelocity.at(1) = x_velocity;

    // Coordinates of the center
    double centerX = center.getX()[0];
    double centerY = center.getX()[1];
    double centerZ = center.getX()[2];

    //full circle
    double maxRadius = radius * distance;

    if (dimension == 2) {
        // 2D Disk
        for (int i = -radius; i <= radius; ++i) {
            for (int j = -radius; j <= radius; ++j) {
                // Calculate the actual position
                double x = centerX + i * distance;
                double y = centerY + j * distance;

                // Calculate displacement vector
                std::array<double, 3> displacement_vector = {x - centerX, y - centerY, 0.0};

                // Calculate the distance from the center using sqrt and pow
                double dist = std::sqrt(std::pow(displacement_vector[0], 2) +
                                        std::pow(displacement_vector[1], 2) +
                                        std::pow(displacement_vector[2], 2));

                // Check if the point is within the circle
                if (dist <= maxRadius) {
                    // Create the particle
                    std::array<double, 3> position = {x, y, centerZ};
                    Particle p(position, addedVelocity, center.getM(),
                        center.getType(), center.getEpsilon(), center.getSigma());

                    // Add the particle to the list
                    allParticles.push_back(p);
                }
            }
        }
    } else if (dimension == 3) {
        // 3D Disk
        for (int i = -radius; i <= radius; ++i) {
            for (int j = -radius; j <= radius; ++j) {
                for (int k = -radius; k <= radius; ++k) {
                    // Calculate the actual position
                    double x = centerX + i * distance;
                    double y = centerY + j * distance;
                    double z = centerZ + k * distance;

                    // Calculate displacement vector
                    std::array<double, 3> displacement_vector = {x - centerX, y - centerY, z - centerZ};

                    // Calculate the distance from the center using sqrt and pow
                    double dist = std::sqrt(std::pow(displacement_vector[0], 2) +
                                            std::pow(displacement_vector[1], 2) +
                                            std::pow(displacement_vector[2], 2));

                    // Check if the point is within the circle
                    if (dist <= maxRadius) {
                        // Create the particle
                        std::array<double, 3> position = {x, y, z};
                        Particle p(position, addedVelocity, center.getM(),
                            center.getType(), center.getEpsilon(), center.getSigma());

                        // Add the particle to the list
                        allParticles.push_back(p);
                    }
                }
            }
        }
    } else {
        // Handle invalid dimension input
        throw std::invalid_argument("Dimension must be 2 or 3");
    }
}


