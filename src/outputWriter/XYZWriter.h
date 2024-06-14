/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "Particle.h"

#include <fstream>
#include <list>

#include "../Container/ParticleContainer.h"

namespace outputWriter {

class XYZWriter {

public:
  XYZWriter();

  virtual ~XYZWriter();

  void plotParticles(BaseParticleContainer &particles, const std::string &filename,
                     int iteration);
};

} // namespace outputWriter
