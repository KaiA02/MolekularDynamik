//
// Created by jh on 04.07.2024.
//

#ifndef MEMBRANEPAIR_H
#define MEMBRANEPAIR_H

#include "../Particle.h"
#include <vector>
/**
 * class MembranePair is used to store the direct and diagonal neighbours of a
 * particle
 */
class MembranePair {
public:
  /**
   * @brief constructor for MembranePair
   * @param center the center particle of the membrane pair
   */
  MembranePair(Particle *center);

  /**
   *
   * @return the center particle of the membrane pair
   */
  Particle *getCenter();

  /**
   * @brief sets the center particle of the membrane pair
   * @param center the center particle of the membrane pair
   */
  void setCenter(Particle *center);

  /**
   * @return the direct neighbours of the center particle
   */
  std::vector<Particle *> getDirect();

  /**
   * @brief adds a particle to the direct neighbours
   * @param p the particle that will be added to the direct neighbours
   */
  void appendDirect(Particle *p);

  /**
   *
   * @return the diagonal neighbours of the center particle
   */
  std::vector<Particle *> getDiagonal();

  /**
   * @brief adds a particle to the diagonal neighbours
   * @param p the particle that will be added to the diagonal neighbours
   */
  void appendDiagonal(Particle *p);

private:
  Particle *center;
  std::vector<Particle *> direct;
  std::vector<Particle *> diagonal;
};
#endif // MEMBRANEPAIR_H
