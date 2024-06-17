#include <gtest/gtest.h>
#include "../Thermostat.h"
#include "../Particle.h"
#include "../Container/LCParticleContainer.h"


TEST(ThermostatTest, Heating) {
    LCParticleContainer pc;
    pc.addParticle(Particle({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1));
    pc.addParticle(Particle({1.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, 2));
    pc.addParticle(Particle({0.0, 1.0, 0.0}, {3.0, 3.0, 3.0}, 3));

    Thermostat thermostat(300.0, 10, 400.0, 10.0);  // Initial temp: 300, target: 400, delta temp: 10

    // Set initial temperature
    thermostat.setInitialTemperature(pc.getParticles());

    double initialTemperature = thermostat.getCurrentTemp(pc.getParticles());
    ASSERT_NEAR(initialTemperature, 300.0, 1.0);  // Initial temperature should be set to 300

    // Apply gradual scaling
    thermostat.gradualScaling(pc.getParticles());

    double currentTemperature = thermostat.getCurrentTemp(pc.getParticles());
    EXPECT_GT(currentTemperature, 300.0);
    EXPECT_LE(currentTemperature, 310.0);  // Max increase is delta_temp
}

TEST(ThermostatTest, Cooling) {
    LCParticleContainer pc;
    pc.addParticle(Particle({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1));
    pc.addParticle(Particle({1.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, 2));
    pc.addParticle(Particle({0.0, 1.0, 0.0}, {3.0, 3.0, 3.0}, 3));

    Thermostat thermostat(500.0, 10, 400.0, 10.0);  // Initial temp: 500, target: 400, delta temp: 10

    // Set initial temperature
    thermostat.setInitialTemperature(pc.getParticles());

    double initialTemperature = thermostat.getCurrentTemp(pc.getParticles());
    ASSERT_NEAR(initialTemperature, 500.0, 1.0);  // Initial temperature should be set to 500

    // Apply gradual scaling
    thermostat.gradualScaling(pc.getParticles());

    double currentTemperature = thermostat.getCurrentTemp(pc.getParticles());
    EXPECT_LT(currentTemperature, 500.0);
    EXPECT_GE(currentTemperature, 490.0);  // Max decrease is delta_temp
}

TEST(ThermostatTest, Holding) {
    LCParticleContainer pc;
    pc.addParticle(Particle({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 1));
    pc.addParticle(Particle({1.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, 2));
    pc.addParticle(Particle({0.0, 1.0, 0.0}, {3.0, 3.0, 3.0}, 3));

    Thermostat thermostat(400.0, 10, 400.0, 10.0);  // Initial temp: 400, target: 400, delta temp: 10

    // Set initial temperature
    thermostat.setInitialTemperature(pc.getParticles());

    double initialTemperature = thermostat.getCurrentTemp(pc.getParticles());
    ASSERT_NEAR(initialTemperature, 400.0, 1.0);  // Initial temperature should be set to 400

    // Apply gradual scaling
    thermostat.gradualScaling(pc.getParticles());

    double currentTemperature = thermostat.getCurrentTemp(pc.getParticles());
    EXPECT_NEAR(currentTemperature, 400.0, 1e-2);  // Should remain near target temperature
}
