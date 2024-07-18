#include "../Particle.h"
#include "../Container/ParticleContainer.h"
#include "spdlog/spdlog.h"
#include <gtest/gtest.h>
#include "../Calculations.h"

TEST(CalculationsTest, TestCalculateLJF) {

ParticleContainer pc;
pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));
pc.addParticle(Particle({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));
pc.addParticle(Particle({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));

// Berechnen der Lennard-Jones-Kräfte
Calculations calculations(pc);
calculations.calculateLJF();

// Überprüfen, ob die Kräfte korrekt berechnet wurden
// Hier können Sie Ihre eigenen Überprüfungen durchführen, basierend auf Ihren Erwartungen
// Ein einfaches Beispiel: (Handcalculatedt values)
std::array<double, 3> f0  = {-120, -120, 0.0};
std::array<double, 3> f1  = {114.375, 5.625, 0.0};
std::array<double, 3> f2  = {5.625, 114.375,0.0};
for (int i = 0; i < 3; ++i) {
EXPECT_NEAR(pc.getParticles().at(0).getF()[i], f0[i], 1e-6);
EXPECT_NEAR(pc.getParticles().at(1).getF()[i], f1[i], 1e-6);
EXPECT_NEAR(pc.getParticles().at(2).getF()[i], f2[i], 1e-6);
}
}

TEST(CalculateSmoothLJFTest, BasicFunctionality) {
// Create particles
Particle* p1 = Particle({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);
Particle* p2 = Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

// Call the function
std::array<double, 3> force = Calculations::calculateSmoothLJF(p1, p2, epsilon, sigma);

// Check the results
EXPECT_NEAR(force[0], -expected_x, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

// Clean up
delete p1;
delete p2;
}

TEST(CalculateSmoothLJFTest, EdgeCaseAtCutoff) {
// Create particles
Particle* p1 = Particle({r_cutoff, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);
Particle* p2 = Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

// Call the function
std::array<double, 3> force = Calculations::calculateSmoothLJF(p1, p2, epsilon, sigma);

// Check the results (expecting zero force)
EXPECT_NEAR(force[0], 0.0, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

// Clean up
delete p1;
delete p2;
}

TEST(CalculateSmoothLJFTest, EdgeCaseAtRL) {
// Create particles
Particle* p1 = Particle({r_l, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);
Particle* p2 = Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

// Call the function
std::array<double, 3> force = Calculations::calculateSmoothLJF(p1, p2, epsilon, sigma);

// Check the results (expecting some force)
EXPECT_NEAR(force[0], -expected_x, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

// Clean up
delete p1;
delete p2;
}

TEST(CalculateSmoothLJFTest, ZeroForceBeyondCutoff) {
// Create particles
Particle* p1 = Particle({r_cutoff + 1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);
Particle* p2 = Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

// Call the function
std::array<double, 3> force = Calculations::calculateSmoothLJF(p1, p2, epsilon, sigma);

// Check the results (expecting zero force)
EXPECT_NEAR(force[0], 0.0, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

// Clean up
delete p1;
delete p2;
}

TEST(CalculateSmoothLJFTest, Correctness) {
// Create particles with known values
Particle* p1 = Particle({1.1, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);
Particle* p2 = Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0);

// Parameters
double epsilon = 1.0;
double sigma = 1.0;

// Known expected result based on manual calculation
double expected_x = -0.00732;

// Call the function
std::array<double, 3> force = Calculations::calculateSmoothLJF(p1, p2, epsilon, sigma);

// Check the results
EXPECT_NEAR(force[0], expected_x, 1e-5);
EXPECT_NEAR(force[1], 0.0, 1e-5);
EXPECT_NEAR(force[2], 0.0, 1e-5);

// Clean up
delete p1;
delete p2;
}



