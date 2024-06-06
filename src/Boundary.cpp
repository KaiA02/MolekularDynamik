//
// Created by joshu on 06.06.2024.
//
#include <array>
#include "Boundary.h"
#include "LinkedCell/Cell.h"

Boundary::Boundary(const std::array<double, 3>& lower, const std::array<double, 3>& upper) : lowerCorner(lower), upperCorner(upper) {}

std::array<double, 3> Boundary::getStartCorner() const {
  return lowerCorner;
}

std::array<double, 3> Boundary::getEndCorner() const {
  return upperCorner;
}

void Boundary::setLowerCorner(const std::array<double, 3>& lower) {
  lowerCorner = lower;
}

void Boundary::setUpperCorner(const std::array<double, 3>& upper) {
  upperCorner = upper;
}

bool Boundary::isCellOnBoundary(Cell c){
        std::array<int, 3> cellId = c.getId();
        if(cellId[0] == lowerCorner[0] || cellId[0] == upperCorner[0] || cellId[1] == lowerCorner[1] || cellId[1] == upperCorner[1] || cellId[2] == lowerCorner[2] || cellId[2] == upperCorner[2]){
            return true;
        }
            return false;
}
