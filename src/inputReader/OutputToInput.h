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

    void readOutput(LCParticleContainer &particleContainer, char *filename);
    bool startsWith(const std::string &str, const std::string &prefix);
    int extractNumberOfPoints(const std::string &line);
};

#endif // OUTPUTTOINPUT_H

