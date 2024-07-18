#ifndef THERMOSTAT_H
#define THERMOSTAT_H
#include "Particle.h"
#include <vector>

/**
 * @brief Thermostat class: used to set/calculate the temperature of the system
 */
class Thermostat {
private:
  double temp_init;
  int n_thermostat;
  double temp_target;
  double delta_temp;

  const static int dimension = 3;

  /**
   * @brief Helper function to calculate the temperature by using the kinetic
   * energy
   * @param particles all particles in the container
   */
  double
  calculateCurrentTemperature(const std::vector<Particle> &particles) const;

  /**
   * @brief Helper function to scalet the velocities by using the scalingFactor
   * scaling factor is calculated by the other methods
   * @param particles all particles in the container
   * @param scalingFactor calculated by the other methods
   */
  void scaleVelocities(std::vector<Particle> &particles,
                       double scalingFactor) const;

public:
  /**
   * @brief constructor for Thermostat
   * @param temp_init initial temperature
   * @param n_thermostat the interval of timesteps when the thermostat is
   * applied
   * @param temp_target target temperature
   * @param delta_temp the difference in temperature that the thermostat can
   * apply
   */
  Thermostat(double temp_init, int n_thermostat, double temp_target,
             double delta_temp);

  /**
   * @return the initial temperature
   */
  const double getTemp_Init() const;

  /**
   * @return the number of timesteps when the thermostat is applied
   */
  const double getN_Thermostat() const;

  /**
   * @return the target temperature
   */
  const double getTemp_Target() const;

  /**
   * @return the difference in temperature that the thermostat can apply
   */
  const double getDelta_Temp() const;

  //@return the current temperature
  const double getCurrentTemp(const std::vector<Particle> &particles) const;

  /**
   * @brief sets the temperature directly
   * calls scaleVelocities()
   * @param particles all particles in the container
   * "setting a temperature directly via velocity scaling"
   */
  void setTemperatureDirectly(std::vector<Particle> &particles) const;

  /**
   * @brief sets the temperature gradualy
   * calls scaleVelocities()
   * @param particles all particles in the container
   * "gradual velocity scaling"
   */
  void gradualScaling(std::vector<Particle> &particles) const;

  /**
   * @brief sets the temperature directly to the initial value
   * calls scaleVelocities()
   * @param particles all particles in the container
   */
  void setInitialTemperature(std::vector<Particle> &particles) const;
};

#endif // THERMOSTAT_H