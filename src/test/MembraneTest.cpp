#include "../Particle.h"
#include "../Container/ParticleContainer.h"
#include "spdlog/spdlog.h"
#include <gtest/gtest.h>
#include "../Calculations.h"
#include "../Membrane/Membrane.h"

TEST(MembraneTest, TestStabilizeMembrane) {
    // Initialize particles
    Particle p1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);
    Particle p2({2.2, 0., 0.0}, {0.0, 0.0, 0.0}, 0, 0);
    Particle p3({0.0, 2.2, 0.0}, {0.0, 0.0, 0.0}, 0, 0);
    Particle p4({2.2, 2.2, 0.0}, {0.0, 0.0, 0.0}, 0, 0); // This will be a diagonal neighbor

    std::vector<Particle*> particles = {&p1, &p2, &p3, &p4};

    // Initialize Membrane with particles, distance and force upwards
    double distance = 2.2;
    double forceUpwards = 0.0;
    Membrane membrane(particles, distance, forceUpwards);

    // Create Calculations instance
    ParticleContainer pc;
    for (auto& p : particles) {
        pc.addParticle(*p);
    }
    Calculations calculations(pc);

    // Define parameters for the harmonic potential
    double stiffness = 10.0; // Example stiffness coefficient
    calculations.setStiffness(stiffness);

    // Stabilize the membrane
    membrane.stabilizeMembrane(calculations);

    // Debug: Output forces for each particle
    std::cout << "Forces after stabilization:" << std::endl;
    std::cout << "p1 force: (" << p1.getF()[0] << ", " << p1.getF()[1] << ", " << p1.getF()[2] << ")" << std::endl;
    std::cout << "p2 force: (" << p2.getF()[0] << ", " << p2.getF()[1] << ", " << p2.getF()[2] << ")" << std::endl;
    std::cout << "p3 force: (" << p3.getF()[0] << ", " << p3.getF()[1] << ", " << p3.getF()[2] << ")" << std::endl;
    std::cout << "p4 force: (" << p4.getF()[0] << ", " << p4.getF()[1] << ", " << p4.getF()[2] << ")" << std::endl;

    // Expected forces based on harmonic potential calculations
    // Manually calculated expected forces
    std::array<double, 3> expectedForceOnP1 = {0.0, 0.0, 0.0};
    std::array<double, 3> expectedForceOnP2 = {0.0, 0.0, 0.0};
    std::array<double, 3> expectedForceOnP3 = {0.0, 0.0, 0.0};
    std::array<double, 3> expectedForceOnP4 = {0.0, 0.0, 0.0};

    // Tolerance for floating-point comparisons
    double tolerance = 1e-6;

    // Check the forces on each particle
    EXPECT_NEAR(p1.getF()[0], expectedForceOnP1[0], tolerance);
    EXPECT_NEAR(p1.getF()[1], expectedForceOnP1[1], tolerance);
    EXPECT_NEAR(p1.getF()[2], expectedForceOnP1[2], tolerance);

    EXPECT_NEAR(p2.getF()[0], expectedForceOnP2[0], tolerance);
    EXPECT_NEAR(p2.getF()[1], expectedForceOnP2[1], tolerance);
    EXPECT_NEAR(p2.getF()[2], expectedForceOnP2[2], tolerance);

    EXPECT_NEAR(p3.getF()[0], expectedForceOnP3[0], tolerance);
    EXPECT_NEAR(p3.getF()[1], expectedForceOnP3[1], tolerance);
    EXPECT_NEAR(p3.getF()[2], expectedForceOnP3[2], tolerance);

    EXPECT_NEAR(p4.getF()[0], expectedForceOnP4[0], tolerance);
    EXPECT_NEAR(p4.getF()[1], expectedForceOnP4[1], tolerance);
    EXPECT_NEAR(p4.getF()[2], expectedForceOnP4[2], tolerance);
}

TEST(MembraneTest, TestApplyMovement) {
    // Initialize particles
    Particle p1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);
    Particle p2({1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}, 0, 0);

    // Create a vector of moving particles
    std::vector<Particle*> movingParticles = {&p1, &p2};

    // Initialize Membrane with particles, distance and force upwards
    double distance = 2.2; // Arbitrary value, not used in applyMovement
    double forceUpwards = 2.0;
    Membrane membrane(movingParticles, distance, forceUpwards);
    membrane.setMovingParticles(movingParticles);


    // Apply movement
    membrane.applyMovement();

    // Expected forces based on the upward force
    std::array<double, 3> expectedForceOnP1 = {0.0, 0.0, 2.0};
    std::array<double, 3> expectedForceOnP2 = {0.0, 0.0, 2.0};

    // Tolerance for floating-point comparisons
    double tolerance = 1e-6;

    // Check the forces on each particle
    EXPECT_NEAR(p1.getF()[0], expectedForceOnP1[0], tolerance);
    EXPECT_NEAR(p1.getF()[1], expectedForceOnP1[1], tolerance);
    EXPECT_NEAR(p1.getF()[2], expectedForceOnP1[2], tolerance);

    EXPECT_NEAR(p2.getF()[0], expectedForceOnP2[0], tolerance);
    EXPECT_NEAR(p2.getF()[1], expectedForceOnP2[1], tolerance);
    EXPECT_NEAR(p2.getF()[2], expectedForceOnP2[2], tolerance);
}