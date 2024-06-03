#include <gtest/gtest.h>
#include "../ParticleContainer.h" // Include your header files
#include "../Particle.h"
#include "../LinkedCell/Cell.h"

TEST(LCParticleContainerTest, GetParticleInNeighbourhood) {
    LCParticleContainer container;
    container.generateCells(3, 3, 3, 1.0); // Generate a 3x3x3 grid

    // Add particles to specific cells
    Particle p1({0.5, 0.5, 0.5}, {0.0, 0.0, 0.0}, 1, 1); // Inside cell (0,0,0)
    Particle p2({1.5, 0.5, 0.5}, {0.0, 0.0, 0.0}, 1, 1); // Inside cell (1,0,0)
    Particle p3({2.5, 0.5, 0.5}, {0.0, 0.0, 0.0}, 1, 1); // Inside cell (2,0,0)
    container.getCellById({0,0,0})->addParticle(p1);
    container.getCellById({1,0,0})->addParticle(p2);
    container.getCellById({2,0,0})->addParticle(p3);

    // Check neighbourhood for cell (1,0,0)
    Cell* c = container.getCellById({1,0,0});
    auto neighbours = container.getParticleInNeighbourhood(*c);

    EXPECT_EQ(neighbours.size(), 3); // p1, p2, p3 should be in the neighbourhood
    EXPECT_NE(std::find(neighbours.begin(), neighbours.end(), p1), neighbours.end());
    EXPECT_NE(std::find(neighbours.begin(), neighbours.end(), p2), neighbours.end());
    EXPECT_NE(std::find(neighbours.begin(), neighbours.end(), p3), neighbours.end());
}

TEST(LCParticleContainerTest, GetCellById) {
    LCParticleContainer container;
    container.generateCells(3, 3, 3, 1.0);

    Cell* cell = container.getCellById({1,1,1});
    std::array<int, 3> res1 = {1,1,1};
    EXPECT_EQ(cell->getId(), res1);

    Cell* nonExistingCell = container.getCellById({10,10,10});
    EXPECT_EQ(nonExistingCell, nullptr);
}

TEST(LCParticleContainerTest, RealocateParticles) {
    LCParticleContainer container;
    container.generateCells(3, 3, 3, 1.0);

    // Add particles
    Particle p1({0.5, 0.5, 0.5}, {0.0, 0.0, 0.0}, 1, 1); // Inside grid
    Particle p2({3.5, 3.5, 3.5}, {0.0, 0.0, 0.0}, 1, 1); // Outside grid
    container.getCellById({0,0,0})->addParticle(p1);
    container.getCellById({2,2,2})->addParticle(p2);

    container.realocateParticles(1);

    // p1 should be reallocated
    Cell* cell = container.getCellById({0,0,0});
    EXPECT_EQ(cell->getParticles().size(), 1);
    std::array<double, 3> res1 = {0.5, 0.5, 0.5};
    EXPECT_EQ(cell->getParticles().at(0).getX(), res1);

    // p2 should be removed
    Cell* cell2 = container.getCellById({2,2,2});
    EXPECT_EQ(cell2->getParticles().size(), 0);
}

TEST(LCParticleContainerTest, GenerateCells) {
    LCParticleContainer container;
    container.generateCells(3, 3, 3, 1.0);

    // Check if cells are generated correctly
    std::array<int, 3> res1 = {0,0,0};
    EXPECT_EQ(container.getCellById({0,0,0})->getId(), res1);
    std::array<int, 3> res2 = {1,1,1};
    EXPECT_EQ(container.getCellById({1,1,1})->getId(), res2);
    std::array<int, 3> res3 = {2,2,2};
    EXPECT_EQ(container.getCellById({2,2,2})->getId(), res3);
    EXPECT_EQ(container.getCellById({3,3,3}), nullptr); // Out of range
}


TEST(LCParticleContainerTest, AddParticle) {
    LCParticleContainer container;
    container.generateCells(3, 3, 3, 1.0); // Generate a 3x3x3 grid

    // Add a particle within the bounds
    Particle p1({0.5, 0.5, 0.5}, {0.0, 0.0, 0.0}, 1, 1); // Inside cell (0,0,0)
    container.addParticle(p1);
    Cell* cell = container.getCellById({0, 0, 0});
    ASSERT_NE(cell, nullptr);
    EXPECT_EQ(cell->getParticles().size(), 1);
    EXPECT_EQ(cell->getParticles().at(0).getX(), p1.getX());

    // Add a particle out of bounds
    Particle p2({3.5, 3.5, 3.5}, {0.0, 0.0, 0.0}, 1, 1); // Outside grid
    container.addParticle(p2);
    Cell* outOfBoundsCell = container.getCellById({3, 3, 3});
    EXPECT_EQ(outOfBoundsCell, nullptr);
    // Check that no cell contains p2
    for (int x = 0; x < 3; ++x) {
        for (int y = 0; y < 3; ++y) {
            for (int z = 0; z < 3; ++z) {
                Cell* c = container.getCellById({x, y, z});
                ASSERT_NE(c, nullptr);
                for (const auto& particle : c->getParticles()) {
                    EXPECT_NE(particle.getX(), p2.getX());
                }
            }
        }
    }
}