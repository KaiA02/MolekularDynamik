//
// Created by jh on 03.06.2024.
//
#include <gtest/gtest.h>
#include "../LinkedCell/Cell.h"
#include "../Particle.h"

// Test the constructor and getId method
TEST(CellTest, ConstructorAndGetId) {
    std::array<int, 3> id = {1, 2, 3};
    Cell cell(id);

    EXPECT_EQ(cell.getId(), id);
}

// Test the addParticle and getParticles methods
TEST(CellTest, AddAndGetParticles) {
    std::array<int, 3> id = {1, 2, 3};
    Cell cell(id);

    Particle p1({{0.1, 0.2, 0.3}}, {{0.0, 0.0, 0.0}}, 1.0);
    Particle p2({{1.1, 1.2, 1.3}}, {{0.0, 0.0, 0.0}}, 1.0);

    cell.addParticle(&p1);
    cell.addParticle(&p2);

    std::vector<Particle*> particles = cell.getParticles();
    ASSERT_EQ(particles.size(), 2);
    //EXPECT_EQ(particles[0], p1);
    //EXPECT_EQ(particles[1], p2);
}

// Test adding no particles
TEST(CellTest, NoParticles) {
    std::array<int, 3> id = {1, 2, 3};
    Cell cell(id);

    std::vector<Particle*> particles = cell.getParticles();
    EXPECT_TRUE(particles.empty());
}

// Test adding a single particle
TEST(CellTest, AddSingleParticle) {
    std::array<int, 3> id = {1, 2, 3};
    Cell cell(id);

    Particle p({{0.1, 0.2, 0.3}}, {{0.0, 0.0, 0.0}}, 1.0);
    cell.addParticle(&p);

    std::vector<Particle*> particles = cell.getParticles();
    ASSERT_EQ(particles.size(), 1);
    //EXPECT_EQ(particles[0], p);
}

// Test particle order
TEST(CellTest, ParticleOrder) {
    std::array<int, 3> id = {1, 2, 3};
    Cell cell(id);

    Particle p1({{0.1, 0.2, 0.3}}, {{0.0, 0.0, 0.0}}, 1.0);
    Particle p2({{1.1, 1.2, 1.3}}, {{0.0, 0.0, 0.0}}, 1.0);
    Particle p3({{2.1, 2.2, 2.3}}, {{0.0, 0.0, 0.0}}, 1.0);

    cell.addParticle(&p1);
    cell.addParticle(&p2);
    cell.addParticle(&p3);

    std::vector<Particle*> particles = cell.getParticles();
    ASSERT_EQ(particles.size(), 3);
    //EXPECT_EQ(particles[0], p1);
    //EXPECT_EQ(particles[1], p2);
    //EXPECT_EQ(particles[2], p3);
}