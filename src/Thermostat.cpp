#include "Thermostat.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <limits>

#include "spdlog/spdlog.h"

Thermostat::Thermostat(double temp_init, int n_thermostat, double temp_target, double delta_temp)
    : temp_init(temp_init), n_thermostat(n_thermostat), temp_target(temp_target), delta_temp(delta_temp) {}

const double Thermostat::getTemp_Init() const {
    return temp_init;
}

const double Thermostat::getN_Thermostat() const {
    return n_thermostat;
}

const double Thermostat::getTemp_Target() const {
    return temp_target;
}

const double Thermostat::getDelta_Temp() const {
    return delta_temp;
}

const double Thermostat::getCurrentTemp(const std::vector<Particle>& particles) const {
    return calculateCurrentTemperature(particles);
}


double Thermostat::calculateCurrentTemperature(const std::vector<Particle>& particles) const {
    double kineticEnergy = 0.0;
    int counter = 0;
    for (const auto& p : particles) {
        std::array<double, 3> comp = {-50.0,-50.0,-50.0};
        if(p.getX() != comp) {
            double skalarproduct = 0;
            for(int i = 0; i < 3; i++) {
                skalarproduct += p.getV().at(i) * p.getV().at(i);
            }

            kineticEnergy += (p.getM() * skalarproduct) / 2.0;
            counter ++;
        }
    }
    spdlog::debug("kinetic Energy: {} ", kineticEnergy);
    return (2.0 * kineticEnergy) / (dimension * counter);

}

void Thermostat::scaleVelocities(std::vector<Particle>& particles, double scalingFactor) const {
    for (auto& p : particles) {
        std::array<double, 3> velocity = {0,0,0};
        for(int i = 0; i < 3; i++) {
            velocity.at(i) += p.getV().at(i) * scalingFactor;
        }
        p.setV(velocity);
    }
}

void Thermostat::setTemperatureDirectly(std::vector<Particle>& particles) const {
    double currentTemperature = calculateCurrentTemperature(particles);
    double scalingFactor = std::sqrt(temp_target / currentTemperature);
    scaleVelocities(particles, scalingFactor);
}

void Thermostat::gradualScaling(std::vector<Particle>& particles) const {
    double currentTemperature = calculateCurrentTemperature(particles);
    double tempChange = std::min(delta_temp, std::abs(temp_target - currentTemperature));
    double temp_new = 0;
    if(temp_target < currentTemperature) {
        temp_new = currentTemperature - tempChange;
    } else {
        temp_new = currentTemperature + tempChange;
    }
    double scalingFactor = std::sqrt(temp_new / currentTemperature);
    scaleVelocities(particles, scalingFactor);
}

void Thermostat::setInitialTemperature(std::vector<Particle> &particles) const {
    double currentTemperature = calculateCurrentTemperature(particles);
    double scalingFactor = std::sqrt(temp_init / currentTemperature);
    scaleVelocities(particles, scalingFactor);
}

