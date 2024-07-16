//
// Created by joshu on 25.05.2024.
//

#ifndef XMLREADER_H
#define XMLREADER_H
#include "../Container/LCParticleContainer.h"
#include "../Container/ParticleContainer.h"
#include "simulation.hxx"

/**
 *@brief class XMLReader is used to read the xml files
 */
class XMLReader {
public:
  XMLReader(const std::string &filePath);
  /**
   *@brief forms the input of the xml files into content
   *@param particleContainer used for the normal ParticleContainer
   */
  void readXML(ParticleContainer &particleContainer);
  /**
   *@brief forms the input of the xml files into content
   *@param particleContainer used for the Linked-Cell ParticleContainer
   */
  void readXML_LC(LCParticleContainer &particleContainer);

  /**
   * @return the type of the particle container
   */
  std::string getParticleContainerType();

  /**
   * @return 0 if we want to have LJF or 1 if we want to have Smooth LJF
   */
  bool getSmoothLJ();

  /**
   * @return the r_L
   */
  double getR_L();

  /**
   * @return check if the thermostat is on or off
   */
  std::string ThermostatON();

  /**
   *
   * @return the initial temperature
   */
  double getTemp_init();

  /**
   *
   * @return the step size of the thermostat
   */
  int getN_Thermostat();

  /**
   *
   * @return the target temperature
   */
  double getTemp_Target();

  /**
   *
   * @return the delta temperature
   */
  double getDelta_Temp();

  /**
   *
   * @return the time step
   */
  std::array<double, 3> getTime();

  /**
   *
   * @return the output type
   */
  std::string getOutputType();

  /**
   *
   * @return the base name
   */
  std::string getBaseName();

  /**
   *
   * @return the write frequency
   */
  int getWriteFrequency();

  /**
   *
   * @return check if the performance measurement is on or off
   */
  xml_schema::boolean getPerformanceMeasurement();

  double getRdfDeltaR();
  /**
   *
   * @return the log level
   */
  std::string getLogLevel();

  /**
   *
   * @return the number of particles in the simulation
   */
  int getNumberOfParticles();

  /**
   *
   * @return the number of cuboids in the simulation
   */
  int getNumberOfCuboids();

  /**
   *
   * @return the number of disks in the simulation
   */
  int getNumberOfDisks();

  xml_schema::boolean getStatisticsOn();

  int getParallelStrategy();

private:
  std::unique_ptr<simulation> sim;
};

#endif // XMLREADER_H
