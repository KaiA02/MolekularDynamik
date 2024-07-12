#include "../LinkedCell/Cell.h"
#include "../Particle.h"
#include "../Container/LCParticleContainer.h"
#include <gtest/gtest.h>
#include <iostream>

TEST(LCParticleContainerTest, GetParticleInNeighbourhood) {
  LCParticleContainer container;
  container.generateCells(3, 3, 3, 1.0); // Generate a 3x3x3 grid

  // Add particles to specific cells
  Particle p1 = Particle ({0.5, 0.5, 0.5}, {0.0, 0.0, 0.0}, 1, 1); // Inside cell (0,0,0)
  Particle p2 = Particle ({1.5, 0.5, 0.5}, {0.0, 0.0, 0.0}, 1, 1); // Inside cell (1,0,0)
  Particle p3 = Particle ({2.5, 0.5, 0.5}, {0.0, 0.0, 0.0}, 1, 1); // Inside cell (2,0,0)
  container.getCellById({0, 0, 0})->addParticle(&p1);
  container.getCellById({1, 0, 0})->addParticle(&p2);
  container.getCellById({2, 0, 0})->addParticle(&p3);

  // Check neighbourhood for cell (1,0,0)
  Cell* c = container.getCellById({1, 0, 0});
  auto neighbours = container.getParticleInNeighbourhood(c->getId());

  EXPECT_EQ(neighbours.size(), 2); // p1, p2, p3 should be in the
  //  neighbourhood
  EXPECT_NE(std::find(neighbours.begin(), neighbours.end(), p1),
            neighbours.end());
  EXPECT_NE(std::find(neighbours.begin(), neighbours.end(), p3),
            neighbours.end());
}

TEST(LCParticleContainerTest, GetCellById) {
  LCParticleContainer container;
  container.generateCells(3, 3, 3, 1.0);

  Cell* cell = container.getCellById({1, 1, 1});
  std::array<int, 3> res1 = {1, 1, 1};
  EXPECT_EQ(cell->getId(), res1);
}

TEST(LCParticleContainerTest, RealocateParticles) {
  LCParticleContainer container;
  container.generateCells(3, 3, 3, 1.0);

  Particle p({3.5, 3.5, 3.5}, {0.0, 0.0, 0.0}, 1, 1); // Outside grid
  container.addParticle(p);
  container.getCellById({2, 2, 2})->addParticle(&p);
  Cell* cell = container.getCellById({2, 2, 2});
  EXPECT_EQ(cell->getParticles().size(), 1);
  container.realocateParticles();
  // p2 should be removed
  cell = container.getCellById({2, 2, 2});
  EXPECT_EQ(cell->getParticles().size(), 0);
}

TEST(LCParticleContainerTest, GenerateCells) {
  LCParticleContainer container;
  container.generateCells(3, 3, 3, 1.0);

  // Check if cells are generated correctly
  std::array<int, 3> res1 = {0, 0, 0};
  EXPECT_EQ(container.getCellById({0, 0, 0})->getId(), res1);
  std::array<int, 3> res2 = {1, 1, 1};
  EXPECT_EQ(container.getCellById({1, 1, 1})->getId(), res2);
  std::array<int, 3> res3 = {2, 2, 2};
  EXPECT_EQ(container.getCellById({2, 2, 2})->getId(), res3);
  // EXPECT_EQ(container.getCellById({3,3,3}), nullptr); // Out of range
}

TEST(LCParticleContainerTest, ApplyGravitation) {
  double g_grav = -12.44;
  LCParticleContainer container;
  container.setG_grav(g_grav);

  std::array<double, 3> position = {1.0, 1.0, 1.0};
  std::array<double, 3> velocity = {0.0, 0.0, 0.0};
  double mass = 1.0;

  // Create and add particles
  Particle p1(position, velocity, mass);
  Particle p2(position, velocity, mass);

  container.addParticle(p1);
  container.addParticle(p2);

  // Apply gravitation
  container.applyGravitation();

  // Check forces
  const auto& particles = container.getParticles();
  ASSERT_EQ(particles.size(), 2);

  for (const auto& p : particles) {
    EXPECT_DOUBLE_EQ(p.getF().at(1), g_grav * mass);
    EXPECT_DOUBLE_EQ(p.getF().at(0), 0.0);
    EXPECT_DOUBLE_EQ(p.getF().at(2), 0.0);
  }
}
