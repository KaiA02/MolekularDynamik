//
// Created by jh on 18.07.2024.
//
#include <gtest/gtest.h>
#include "../Membrane/Membrane.h"

TEST(MembraneTest, applyMovement) {
    Particle p1, p2, p3, p4;
    std::vector<Particle> real_particles = {p1, p2, p3, p4};
    std::vector<Particle*> particles = {&real_particles.at(0), &real_particles.at(1), &real_particles.at(2), &real_particles.at(3)};
    Membrane membrane(particles, 2.2, 10.0);
    membrane.setMovingParticles(particles);

    // Initialize the membrane nodes with some values
    membrane.applyMovement();


    // Verify the force application logic
    for (auto p : real_particles) {
        EXPECT_DOUBLE_EQ(p.getF()[0], 0.0);
        EXPECT_DOUBLE_EQ(p.getF()[1], 0.0);
        EXPECT_DOUBLE_EQ(p.getF()[2], 10.0);
    }
}
