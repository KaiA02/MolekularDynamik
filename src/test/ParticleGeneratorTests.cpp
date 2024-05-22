#include "../Particle.h"
#include "spdlog/spdlog.h"
#include <gtest/gtest.h>
#include "../ParticleGenerator.h"
#include "../utils/MaxwellBoltzmannDistribution.h"

TEST(ParticleGeneratorTest, GenerateCuboidTest) {
    ParticleGenerator generator;
    Particle start({0, 0, 0}, {0, 0, 0}, 1.0);
    int n1 = 40, n2 = 8, n3 = 1;
    double distance = 1.1225;
    double meanVelocity = 0.1;
   // std::array<double, 3> maxwellVelocity = maxwellBoltzmannDistributedVelocity(meanVelocity, 3);

    generator.generateCuboid(start, n1, n2, n3, distance, meanVelocity,3);
    std::vector<Particle> cubeParticles = generator.getCube();

    // Überprüfen, ob die Anzahl der generierten Partikel korrekt ist
    ASSERT_EQ(cubeParticles.size(), n1 * n2 * n3);

    // Überprüfen, ob die Partikel korrekt generiert wurden
    //for (int x = 0; x < n1; ++x) {
    // for (int y = 0; y < n2; ++y) {
    //    for (int z = 0; z < n3; ++z) {
    //       Particle expectedParticle(
    //              {start.getX()[0] + x * distance, start.getX()[1] + y * distance, start.getX()[2] + z * distance},
    //             {start.getV()[0] + maxwellVelocity[0], start.getV()[1] + maxwellVelocity[1], start.getV()[2] + maxwellVelocity[2]},
    //            start.getM(), start.getType());

    //     int index = x * n2 * n3 + y * n3 + z;
    //   EXPECT_EQ(cubeParticles[index], expectedParticle);
    //  }
    //  }
    // }
}
