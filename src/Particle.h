/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>

class Particle {

private:
  /**
   * Position of the particle
   */
  std::array<double, 3> x;

  /**
   * Velocity of the particle
   */
  std::array<double, 3> v;

  /**
   * Force effective on this particle
   */
  std::array<double, 3> f;

  /**
   * Force which was effective on this particle
   */
  std::array<double, 3> old_f;

  /**
   * Mass of this particle
   */
  double m;

  /**
   * Type of the particle. Use it for whatever you want (e.g. to separate
   * molecules belonging to different bodies, matters, and so on)
   */
  int type;

  /**
   * Lennard-Jones parameter epsilon
   */
  double epsilon;

  /**
   * Lennard-Jones parameter sigma
   */
  double sigma;

  bool isHalo;
  std::array<int, 3> ID;

public:
  explicit Particle(int type = 0);

  Particle(const Particle &other);

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
      int type = 0, double epsilon = 5, double sigma = 1);

  ~Particle();

  // virtual ~Particle();

  const std::array<double, 3> &getX() const;

  const std::array<double, 3> &getV() const;

  const std::array<double, 3> &getF() const;

  const std::array<double, 3> &getOldF() const;

  double getM() const;

  int getType() const;

  /**
   *
   * @return the Lennard-Jones parameter epsilon
   */
  double getSigma() const;

  /**
   *
   * @return the Lennard-Jones parameter sigma
   */
  double getEpsilon() const;

  bool operator==(Particle &other);

  std::string toString() const;

  void setX(const std::array<double, 3> &newPosition);

  void setV(const std::array<double, 3> &newVelocity);

  void setF(const std::array<double, 3> &newForce);

  void setM(const double &newMass);

  void setType(const int &newType);

  /**
   * @brief sets the Lennard-Jones parameter sigma
   * @param sigma the Lennard-Jones parameter sigma
   */
  void setSigma(const double &sigma);

  /**
   * sets the Lennard-Jones parameter epsilon
   * @param epsilon the Lennard-Jones parameter epsilon
   */
  void setEpsilon(const double &epsilon);

  /**
   * @brief sets the particle to a halo particle
   * @param halo true if the particle is a halo particle, false otherwise
   */
  void setIsHalo(bool halo);

  /**
   * @return true if the particle is a halo particle, false otherwise
   */
  bool getIsHalo();

  /**
   *
   * @return the id of the particle
   */
  std::array<int, 3> getID();

  /**
   * @brief sets the id of the particle
   * @param id the id of the particle
   */
  void setID(std::array<int, 3> id);

  bool operator==(const Particle &other) const;

  /**
   * @brief parks the particle at position -50 -50 -50 with velocity and force 0
   */
  void park();

  /**
   * @brief checks if the particle n is a neighbour of this particle
   * @param n the particle to check if it is a neighbour
   * @return true if n is a neighbour of this particle, false otherwise
   */
  bool isNeighbour(Particle *n);
};

std::ostream &operator<<(std::ostream &stream, Particle &p);
