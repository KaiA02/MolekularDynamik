//
// Created by joshu on 25.05.2024.
//

#ifndef XMLREADER_H
#define XMLREADER_H

#include "../Container/ParticleContainer.h"
#include "simulation.hxx"
#include <iostream>

#include "Container/LCParticleContainer.h"

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
  std::string getParticleContainerType();
  std::string getInputType();
  std::array<double, 3> getTime();
  std::string getOutputType();
  std::string getBaseName();
  int getWriteFrequency();
  xml_schema::boolean getPerformanceMeasurement();
  std::string getLogLevel();
  int getNumberOfParticles();
  int getNumberOfCuboids();
  int getNumberOfDisks();

private:
  std::unique_ptr<simulation> sim;
};

#endif // XMLREADER_H
