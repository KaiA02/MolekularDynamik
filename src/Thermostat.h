#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#include <vector>

#include "Particle.h"

class Thermostat {
private:
    double temp_init;
    int n_thermostat;
    double temp_target;
    double delta_temp;

    int dimension;

    // Helper function to calculate current temperature
    double calculateCurrentTemperature(const std::vector<Particle>& particles) const;

    // Helper function to scale velocities
    void scaleVelocities(std::vector<Particle>& particles, double scalingFactor) const;

public:
    Thermostat(double temp_init, int n_thermostat, double temp_target, double delta_temp, int dimension);

    const double getTemp_Init() const;
    const double getN_Thermostat() const;
    const double getTemp_Target() const;
    const double getDelta_Temp() const;
    const double getCurrentTemp(const std::vector<Particle>& particles) const;

    // Apply thermostat periodically (uses graudal velocity scaling)
    //First option
    void applyThermostatPeriodically(std::vector<Particle>& particles, int currentStep) const;

    // Set temperature directly via velocity scaling
    //Second Option
    void setTemperatureDirectly(std::vector<Particle>& particles) const;

    // Gradual velocity scaling
    void gradualScaling(std::vector<Particle>& particles) const;
};

#endif //THERMOSTAT_H