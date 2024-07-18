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
pc.addParticle(Particle({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));

// Create Calculations instance
Calculations calculations(pc);

// Define epsilon and sigma for the Lennard-Jones potential
double epsilon = 1.0; // Example value
double sigma = 1.0;   // Example value

// Define the cutoff radius
double r_cutoff = 2.5 * sigma;

// Assign the cutoff radius to the Calculations instance if necessary
calculations.setCutoffRadius(r_cutoff); // Assuming such a method exists

// Calculate Lennard-Jones forces for all particle pairs
for (auto& pc.getParticles()[0] : pc.getParticles()) {
for (auto& pc.getParticles()[1] : pc.getParticles()) {
if (&pc.getParticles()[0] != &pc.getParticles()[1]) {
std::array<double, 3> force = calculations.calculateLJF(&pc.getParticles()[0], &pc.getParticles()[1], epsilon, sigma);
std::array<double, 3> currentForce = pc.getParticles()[0].getF();
pc.getParticles()[0].setF({currentForce[0] + force[0], currentForce[1] + force[1], currentForce[2] + force[2]});
}
}
}

// Expected forces (precomputed or manually calculated)
std::array<double, 3> f0 = {-24.0, -24.0, 0.0}; // Adjusted to realistic values
std::array<double, 3> f1 = {12.0, 12.0, 0.0};   // Adjusted to realistic values
std::array<double, 3> f2 = {12.0, 12.0, 0.0};   // Adjusted to realistic values

// Tolerance for floating-point comparisons
double tolerance = 1e-6;

// Check the forces on each particle
for (int i = 0; i < 3; ++i) {
EXPECT_NEAR(pc.getParticles().at(0).getF()[i], f0[i], tolerance);
EXPECT_NEAR(pc.getParticles().at(1).getF()[i], f1[i], tolerance);
EXPECT_NEAR(pc.getParticles().at(2).getF()[i], f2[i], tolerance);
}
}

TEST(CalculateSmoothLJFTest, BasicFunctionality) {
// Create particles
ParticleContainer pc;
pc.addParticle(Particle({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));
pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

Calculations calculations(pc);

// Call the function
std::array<double, 3> force = calculations.calculateSmoothLJF(pc.getParticles()[0], pc.getParticles()[1], epsilon, sigma);

// Check the results
EXPECT_NEAR(force[0], -expected_x, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

}

TEST(CalculateSmoothLJFTest, EdgeCaseAtCutoff) {
// Create particles
ParticleContainer pc;
pc.addParticle(Particle({2.6, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0)); //r_cutoff=2.6
pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

Calculations calculations(pc);

// Call the function
std::array<double, 3> force = calculations.calculateSmoothLJF(pc.getParticles()[0], pc.getParticles()[1], epsilon, sigma);

// Check the results (expecting zero force)
EXPECT_NEAR(force[0], 0.0, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

}

TEST(CalculateSmoothLJFTest, EdgeCaseAtRL) {
ParticleContainer pc;
pc.addParticle(Particle({r_l, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));
pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

Calculations calculations(pc);

// Call the function
std::array<double, 3> force = calculations.calculateSmoothLJF(pc.getParticles()[0], pc.getParticles()[1], epsilon, sigma);

// Check the results (expecting some force)
EXPECT_NEAR(force[0], -expected_x, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

}

TEST(CalculateSmoothLJFTest, ZeroForceBeyondCutoff) {
// Create particles
ParticleContainer pc;
pc.addParticle(Particle({3.6, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0)); //r_cutoff= 2.6 + 1
pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

// Call the function
Calculations calc(pc);
std::array<double, 3> force = calc.calculateSmoothLJF(pc.getParticles()[0], pc.getParticles()[1], epsilon, sigma);

// Check the results (expecting zero force)
EXPECT_NEAR(force[0], 0.0, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

}

TEST(CalculateSmoothLJFTest, Correctness) {
// Create particles with known values
ParticleContainer pc;
pc.addParticle(Particle({1.1, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0)();
pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

// Known expected result based on manual calculation
double expected_x = -0.00732;

// Call the function
Calculations calc(pc);
std::array<double, 3> force = calc.calculateSmoothLJF(pc.getParticles()[0], pc.getParticles()[1], epsilon, sigma);

// Check the results
EXPECT_NEAR(force[0], expected_x, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

}


