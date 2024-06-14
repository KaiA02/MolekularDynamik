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

    generator.generateCuboid(start, n1, n2, n3, distance, meanVelocity, 3, 40);
    std::vector<Particle> cubeParticles = generator.getAllParticles();

    // Check whether the number of particles generated is correct
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


TEST(ParticleGeneratorTest, Generate2DDiskTest) {
    Particle center({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0.0);
    ParticleGenerator generator;
    int radius = 2;
    double distance = 1.0;
    int dimension = 2;
    double temp_init = 40;
    double meanVelocity = 0.1;

    generator.generateDisk(center, radius, distance, meanVelocity, dimension, temp_init);
    std::vector<Particle> disk = generator.getAllParticles();

    // The expected number of particles in a 2D disk with radius 2 is 13 (1 + 4 + 8)
    ASSERT_EQ(disk.size(), 13);

    // Verify the positions of the particles
    for (const auto& particle : disk) {
        double x = particle.getX()[0];
        double y = particle.getX()[1];
        double z = particle.getX()[2];

        double dist = std::sqrt(x * x + y * y);
        EXPECT_LE(dist, 2.0);
        EXPECT_EQ(z, 0.0);
    }
}

TEST(ParticleGeneratorTest, Generate3DDiskTest) {
    Particle center({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0.0);
    ParticleGenerator generator;
    int radius = 1;
    double distance = 1.0;
    int dimension = 3;
    double temp_init = 40;
    double meanVelocity = 0.1;

    generator.generateDisk(center, radius, distance, meanVelocity, dimension, temp_init);
    std::vector<Particle> disk = generator.getAllParticles();

    // The expected number of particles in a 3D disk with radius 1 is 7 (1 + 6)
    ASSERT_EQ(disk.size(), 7);

    // Verify the positions of the particles
    for (const auto& particle : disk) {
        double x = particle.getX()[0];
        double y = particle.getX()[1];
        double z = particle.getX()[2];

        double dist = std::sqrt(x * x + y * y + z * z);
        EXPECT_LE(dist, 1.0);
    }
}

TEST(ParticleGeneratorTest, InvalidDimension_Disk) {
    Particle center;
    ParticleGenerator generator;

    // Expect an invalid_argument exception for invalid dimension
    EXPECT_THROW(generator.generateDisk(center, 2, 1.0, 0.1, 4, 0), std::invalid_argument);
}
