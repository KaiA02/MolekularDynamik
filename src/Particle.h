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
   * Lennard-Jones parameter sigma
   */
 double sigma;

 /**
   * Lennard-Jones parameter epsilon
   */
 double epsilon;


 bool isHalo;

public:
  explicit Particle(int type = 0);

  Particle(const Particle &other);

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
      int type = 0, double epsilon = 5, double sigma = 1);

  ~Particle();

  //virtual ~Particle();

  const std::array<double, 3> &getX() const;

  const std::array<double, 3> &getV() const;

  const std::array<double, 3> &getF() const;

  const std::array<double, 3> &getOldF() const;

  double getM() const;

  int getType() const;

  double getSigma() const;
  double getEpsilon() const;

  bool operator==(Particle &other);

  std::string toString() const;

  void setX(const std::array<double, 3> &newPosition);

  void setV(const std::array<double, 3> &newVelocity);

  void setF(const std::array<double, 3> &newForce);

 void setM(const double &newMass);

 void setType(const int &newType);

 void setSigma(const double &sigma);
 void setEpsilon(const double &epsilon);

 void setIsHalo(bool halo);

 bool getIsHalo();

 bool operator==(const Particle &other) const;
 void park();

};

std::ostream &operator<<(std::ostream &stream, Particle &p);
