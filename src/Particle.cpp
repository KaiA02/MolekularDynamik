/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"
#include <iostream>

Particle::Particle(int type_arg) {
  type = type_arg;
  //std::cout << "Particle generated!" << std::endl;
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
}

Particle::Particle(const Particle &other) {
  x = other.x;
  v = other.v;
  f = other.f;
  old_f = other.old_f;
  m = other.m;
  type = other.type;
  sigma = other.sigma;
  epsilon = other.epsilon;
  spdlog::debug("Particle generated by copy");
}

// Todo: maybe use initializater list instead of copy?
Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, int type_arg, double e, double s) {
  x = x_arg;
  v = v_arg;
  m = m_arg;
  type = type_arg;
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
  epsilon = e;
  sigma = s;
  spdlog::debug("Particle generated");
}

Particle::~Particle() {
  spdlog::debug("Particle destroyed");
}


const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

const std::array<double, 3> &Particle::getF() const { return f; }

const std::array<double, 3> &Particle::getOldF() const { return old_f; }

double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

double Particle::getSigma() const { return sigma; }
double Particle::getEpsilon() const { return epsilon; }

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << f
         << " old_f: " << old_f << " type: " << type;
  return stream.str();
}

bool Particle::operator==(Particle &other) {
  return (x == other.x) and (v == other.v) and (f == other.f) and
         (type == other.type) and (m == other.m) and (old_f == other.old_f);
}

std::ostream &operator<<(std::ostream &stream, Particle &p) {
  stream << p.toString();
  return stream;
}

void Particle::setX(const std::array<double, 3> &newPosition) {
  x = newPosition;
}

void Particle::setV(const std::array<double, 3> &newVelocity) {
  v = newVelocity;
}

void Particle::setF(const std::array<double, 3> &newForce) {
  old_f = f;
  f = newForce;
}

void Particle::setM(const double &newMass) {
  m = newMass;
}
void Particle::setType(const int &newType) {
  type = newType;
}

void Particle::setSigma(const double &s) {
  sigma = s;
}
void Particle::setEpsilon(const double &e) {
  epsilon = e;
}

bool Particle::operator==(const Particle &other) const {
  // Vergleiche die Attribute der beiden Partikel
  if (this->getX() != other.getX()) return false;
  if (this->getV() != other.getV()) return false;
  if (this->getF() != other.getF()) return false;
  if (this->getOldF() != other.getOldF()) return false;
  if (this->getM() != other.getM()) return false;
  if (this->getType() != other.getType()) return false;
  if (this->getSigma() != other.getSigma()) return false;
  if (this->getEpsilon() != other.getEpsilon()) return false;

  // Wenn alle Attribute übereinstimmen, gib true zurück
  return true;
}

void Particle::park() {
  setF({0.0,0.0,0.0});
  setX({-50.0,-50.0,-50.0});
  setV({0.0,0.0,0.0});
  setM(0);
  setType(0);
}

bool Particle::getIsHalo() {
  return isHalo;
}

 void Particle::setIsHalo(bool halo) {
   isHalo = halo;
 }





