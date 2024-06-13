//
// Created by kaiarenja on 13.06.24.
//

#ifndef BASEPARTICLECONTAINER_H
#define BASEPARTICLECONTAINER_H
#include <vector>


class Cell;
class Particle;
/**
 * @brief The BaseParticleContainer class
 * This class is a base container for Particles. It contains a vector of Particles
 * and provides basic operations for managing Particles.
 */
class BaseParticleContainer {
public:
    virtual ~BaseParticleContainer() = default;
    /**
    * @brief add a Particle to the vector
    * @param particle Particle to add
    */
    virtual void addParticle(const Particle &particle) = 0;
    virtual std::vector<Cell> getCells() = 0;
    /**
     * @brief replaces a particle in the vector particles at index position
     * @param p new Particle
     * @param position index in particles
     */
    virtual void setParticle(Particle p, int position) = 0;

    /**
     * @brief resets all particles in the container
     */
    virtual void resetParticles() = 0;

    /**
     * @brief gets all particles in the container
     * @return reference to vector of particles
     */
    virtual std::vector<Particle>& getParticles() = 0;
    virtual const std::vector<Particle> &getParticles() const = 0;

    /**
     * @brief gets the number of particles in the container
     * @return number of particles
     */
    virtual int size() const = 0;

    virtual std::vector<Particle>::iterator begin() = 0;
    virtual std::vector<Particle>::iterator end() = 0;
    virtual std::vector<Particle>::const_iterator begin() const = 0;
    virtual std::vector<Particle>::const_iterator end() const = 0;
    virtual void handleLJFCalculation() = 0;

protected:
    std::vector<Particle> particles;
};



#endif //BASEPARTICLECONTAINER_H
