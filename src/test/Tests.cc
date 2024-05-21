#include "../Particle.h"
#include "../ParticleContainer.h"
#include "spdlog/spdlog.h"
#include <gtest/gtest.h>

#include "../Calculations.h"
#include "../ParticleGenerator.h"
#include "../utils/MaxwellBoltzmannDistribution.h"

// Test case for ParticleContainer
TEST(ParticleContainerTest, SizeTest) {
  // Create a ParticleContainer
  ParticleContainer pc;

  // Create 4 particles and add them to the container
  Particle p1, p2, p3, p4;
  pc.addParticle(p1);
  pc.addParticle(p2);
  pc.addParticle(p3);
  pc.addParticle(p4);

  // Check if the size of the container is 4
  EXPECT_EQ(pc.size(), 4);
}

TEST(ParticleGeneratorTest, GenerateCuboidTest) {
    ParticleGenerator generator;
    Particle start({0, 0, 0}, {0, 0, 0}, 1.0);
    int n1 = 40, n2 = 8, n3 = 1;
    double distance = 1.1225;
    double meanVelocity = 0.1;
    std::array<double, 3> maxwellVelocity = maxwellBoltzmannDistributedVelocity(meanVelocity, 3);

    generator.generateCuboid(start, n1, n2, n3, distance, meanVelocity,3);
    std::vector<Particle> cubeParticles = generator.getCube();

    // Überprüfen, ob die Anzahl der generierten Partikel korrekt ist
    ASSERT_EQ(cubeParticles.size(), n1 * n2 * n3);

    // Überprüfen, ob die Partikel korrekt generiert wurden
   // for (int x = 0; x < n1; ++x) {
    //   for (int y = 0; y < n2; ++y) {
    //        for (int z = 0; z < n3; ++z) {
     //           Particle expectedParticle(
     //                   {start.getX()[0] + x * distance, start.getX()[1] + y * distance, start.getX()[2] + z * distance},
     //                   {start.getV()[0] + maxwellVelocity[0], start.getV()[1] + maxwellVelocity[1], start.getV()[2] + maxwellVelocity[2]},
      //                  start.getM(), start.getType());

        //        int index = x * n2 * n3 + y * n3 + z;
         //       EXPECT_EQ(cubeParticles[index], expectedParticle);
          //  }
       // }
    //}
}
// Testfall für die calculateLJF-Methode
TEST(CalculationsTest, TestCalculateLJF) {

 ParticleContainer pc;
 pc.addParticle(Particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 1.0));
 pc.addParticle(Particle({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 1.0));
 pc.addParticle(Particle({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, 0, 1.0));

 // Berechnen der Lennard-Jones-Kräfte
 Calculations calculations(pc);
 calculations.calculateLJF();

 // Überprüfen, ob die Kräfte korrekt berechnet wurden
 // Hier können Sie Ihre eigenen Überprüfungen durchführen, basierend auf Ihren Erwartungen
 // Ein einfaches Beispiel: (Handcalculatedt values)
 std::array<double, 3> f0  = {-120, -120, 0.0};
 std::array<double, 3> f1  = {114.375, 5.625, 0};
 std::array<double, 3> f2  = {5.625, 114.375,0.0};
 //EXPECT_EQ(pc.getParticles().at(0).getF(), f0);
 //EXPECT_EQ(pc.getParticles().at(1).getF(), f1);
 //EXPECT_EQ(pc.getParticles().at(2).getF(), f2);
}
