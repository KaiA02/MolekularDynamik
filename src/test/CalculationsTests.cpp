#include "../Particle.h"
#include "../Container/ParticleContainer.h"
#include "spdlog/spdlog.h"
#include <gtest/gtest.h>
#include "../Calculations.h"

TEST(CalculationsTest, TestCalculateLJF) {
    // Initialize ParticleContainer and add particles
    ParticleContainer pc;
    pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));
    pc.addParticle(Particle({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));

    // Create Calculations instance
    Calculations calculations(pc);

    // Define epsilon and sigma for the Lennard-Jones potential
    double epsilon = 1.0; // Example value
    double sigma = 1.0;   // Example value

    // Define the cutoff radius
    double r_cutoff = 2.5 * sigma;

    // Assign the cutoff radius to the Calculations instance if necessary
    calculations.setR_cutoff(r_cutoff);

    // Calculate Lennard-Jones forces for all particle pairs
    for (auto& p1 : pc.getParticles()) {
        for (auto& p2 : pc.getParticles()) {
            if (&p1 != &p2) {
                std::array<double, 3> force = calculations.calculateLJF(&p1, &p2, epsilon, sigma);
                std::array<double, 3> currentForce = p1.getF();
                p1.setF({currentForce[0] + force[0], currentForce[1] + force[1], currentForce[2] + force[2]});
            }
        }
    }

    // Expected forces
    std::array<double, 3> f0 = {-24.0, 0.0, 0.0}; // Adjusted to realistic values
    std::array<double, 3> f1 = {24.0, 0.0, 0.0};  // Adjusted to realistic values

    // Tolerance for floating-point comparisons
    double tolerance = 1e-6;

    // Check the forces on each particle
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(pc.getParticles().at(0).getF()[i], f0[i], tolerance);
        EXPECT_NEAR(pc.getParticles().at(1).getF()[i], f1[i], tolerance);
    }
}

TEST(CalculationsTest, TestcalculatsmoothLJF) {
    // Initialize ParticleContainer and add particles
    ParticleContainer pc;
    pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));
    pc.addParticle(Particle({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));

    // Create Calculations instance
    Calculations calculations(pc);

    // Define epsilon and sigma for the Lennard-Jones potential
    double epsilon = 1.0; // Example value
    double sigma = 1.0;   // Example value

    // Define the cutoff radius
    double r_cutoff = 2.5 * sigma;

    // Assign the cutoff radius to the Calculations instance if necessary
    calculations.setR_cutoff(r_cutoff);

    // Calculate Lennard-Jones forces for all particle pairs
    for (auto& p1 : pc.getParticles()) {
        for (auto& p2 : pc.getParticles()) {
            if (&p1 != &p2) {
                std::array<double, 3> force = calculations.calculateSmoothLJF(&p1, &p2, epsilon, sigma);
                std::array<double, 3> currentForce = p1.getF();
                p1.setF({currentForce[0] + force[0], currentForce[1] + force[1], currentForce[2] + force[2]});
            }
        }
    }

    // Expected forces
    std::array<double, 3> f0 = {-15.552, 0.0, 0.0}; // Adjusted to realistic values
    std::array<double, 3> f1 = {15.552, 0.0, 0.0};  // Adjusted to realistic values

    // Tolerance for floating-point comparisons
    double tolerance = 1e-6;

    // Check the forces on each particle
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(pc.getParticles().at(0).getF()[i], f0[i], tolerance);
        EXPECT_NEAR(pc.getParticles().at(1).getF()[i], f1[i], tolerance);
    }
}



