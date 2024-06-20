//
// Created by jh on 18.06.2024.
//

#ifndef OUTPUTTOINPUT_H
#define OUTPUTTOINPUT_H

#include "../Container/LCParticleContainer.h"
#include <string>

#include "../Particle.h"

#include <list>

#include "../Container/ParticleContainer.h"
/**
 * @brief FileReader class
 */
class OutputToInput {

public:
  OutputToInput();
  virtual ~OutputToInput();

  /**
   * reads the output file and adds the resulting particles to the
   * particleContainer
   * @param particleContainer the container to which the particles will be
   * added
   * @param filename the name of the output file
   */
  void readOutput(LCParticleContainer &particleContainer, char *filename);

  /**
   * @brief checks if a string starts with a prefix
   * @param str the string to check
   * @param prefix the prefix to check for
   * @return true if the string starts with the prefix, false otherwise
   */
  bool startsWith(const std::string &str, const std::string &prefix);
  int extractNumberOfPoints(const std::string &line);
};

#endif // OUTPUTTOINPUT_H
