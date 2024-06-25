#include "../Particle.h"
#include "../Container/ParticleContainer.h"
#include "spdlog/spdlog.h"
#include <gtest/gtest.h>
#include "../Calculations.h"

TEST(CalculationsTest, TestCalculateLJF) {

    ParticleContainer pc;
    pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0));
    pc.addParticle(Particle({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0));
    pc.addParticle(Particle({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, 0));

    // Berechnen der Lennard-Jones-Kräfte
    Calculations calculations(pc);
    calculations.calculateLJF();

    // Überprüfen, ob die Kräfte korrekt berechnet wurden
    // Hier können Sie Ihre eigenen Überprüfungen durchführen, basierend auf Ihren Erwartungen
    // Ein einfaches Beispiel: (Handcalculatedt values)
    std::array<double, 3> f0  = {-120, -120, 0.0};
    std::array<double, 3> f1  = {114.375, 5.625, 0.0};
    std::array<double, 3> f2  = {5.625, 114.375,0.0};
    std::array<double, 3> f00 = pc.getParticles().at(0).getF();
    std::array<double, 3> f01 = pc.getParticles().at(1).getF();
    std::array<double, 3> f02 = pc.getParticles().at(2).getF();
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(f00[i], f0[i], 1e-6);
        EXPECT_NEAR(f01[i], f1[i], 1e-6);
        EXPECT_NEAR(f02[i], f2[i], 1e-6);
    }
}
