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
                                       double meanVelocity, int dimension) {
  std::array<double, 3> maxwellVelocity =
      maxwellBoltzmannDistributedVelocity(meanVelocity, dimension);

  for (int x = 0; x < n1; x++) {
    for (int y = 0; y < n2; y++) {
      for (int z = 0; z < n3; z++) {
        std::array<double, 3> addedVelocity;
        for (int i = 0; i < 3; i++) {
          addedVelocity[i] = start.getV()[i] + maxwellVelocity[i];
        }

        Particle p({start.getX()[0] + x * distance,
                    start.getX()[1] + y * distance,
                    start.getX()[2] + z * distance},
                   addedVelocity, start.getM(), start.getType());
        cube.push_back(p);
      }
    }
  }
}

std::vector<Particle> ParticleGenerator::getCube() { return cube; }

void ParticleGenerator::generateDisk(const Particle &center, int radius, double distance, int dimension) {
  // Clear previous particles
    disk.clear();

    // Coordinates of the center
    double centerX = center.getX()[0];
    double centerY = center.getX()[1];
    double centerZ = center.getX()[2];

    // Initial velocity of the center
    std::array<double, 3> initialVelocity = center.getV();

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
                if (dist <= radius * distance) {
                    // Create the particle
                    std::array<double, 3> position = {x, y, centerZ};
                    Particle p(position, initialVelocity, center.getM(), center.getType());

                    // Add the particle to the list
                    disk.push_back(p);
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
                    if (dist <= radius * distance) {
                        // Create the particle
                        std::array<double, 3> position = {x, y, z};
                        Particle p(position, initialVelocity, center.getM(), center.getType());

                        // Add the particle to the list
                        disk.push_back(p);
                    }
                }
            }
        }
    } else {
        // Handle invalid dimension input
        throw std::invalid_argument("Dimension must be 2 or 3");
    }
}

std::vector<Particle> ParticleGenerator::getDisk() { return disk;}
