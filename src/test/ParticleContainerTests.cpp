#include "../Particle.h"
#include "../Container/ParticleContainer.h"
#include "spdlog/spdlog.h"
#include <gtest/gtest.h>

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


