//
// Created by joshu on 06.06.2024.
//
#include "LinkedCell/Cell.h"
#include <array>

#ifndef BOUNDARY_H
#define BOUNDARY_H
class Boundary {
private:
  std::array<double, 3> lowerCorner;
  std::array<double, 3> upperCorner;

public:
  // Constructor
  Boundary(const std::array<double, 3> &lower,
           const std::array<double, 3> &upper);

  // Getters
  std::array<double, 3> getStartCorner() const;
  std::array<double, 3> getEndCorner() const;

  // Setters
  void setLowerCorner(const std::array<double, 3> &lower);
  void setUpperCorner(const std::array<double, 3> &upper);

  bool isCellOnBoundary(Cell c);
};
#endif // BOUNDARY_H
