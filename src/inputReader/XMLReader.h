//
// Created by joshu on 25.05.2024.
//

#ifndef XMLREADER_H
#define XMLREADER_H

#include "../ParticleContainer.h"
#include "simulation.hxx"
#include <iostream>

class XMLReader {
public:
  XMLReader(const std::string &filePath);
  void readXML(ParticleContainer &particleContainer);
  std::string getInputType();
  std::array<double, 3> getTime();
  std::string getOutputType();
  std::string getBaseName();
  int getWriteFrequency();
  xml_schema::boolean getPerformanceMeasurement();
  std::string getLogLevel();
  int getNumberOfParticles();
  int getNumberOfCuboids();

private:
  std::unique_ptr<simulation> sim;
};

#endif // XMLREADER_H
