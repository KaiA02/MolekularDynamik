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

private:
  std::unique_ptr<simulation> sim;
};

#endif // XMLREADER_H
