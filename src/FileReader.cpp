/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "ParticleGenerator.h"

FileReader::FileReader() = default;

FileReader::~FileReader() = default;

void FileReader::readFile(ParticleContainer& particleContainer, char *filename) {
  std::array<double, 3> x;
  std::array<double, 3> v;
  double m;
  int num_particles = 0;

  std::ifstream input_file(filename);
  std::string tmp_string;

  if (input_file.is_open()) {

    getline(input_file, tmp_string);
    std::cout << "Read line: " << tmp_string << std::endl;

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << std::endl;
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    std::cout << "Reading " << num_particles << "." << std::endl;
    getline(input_file, tmp_string);
    std::cout << "Read line: " << tmp_string << std::endl;

    for (int i = 0; i < num_particles; i++) {
      std::istringstream datastream(tmp_string);

      for (auto &xj : x) {
        datastream >> xj;
      }
      for (auto &vj : v) {
        datastream >> vj;
      }
      if (datastream.eof()) {
        std::cout
            << "Error reading file: eof reached unexpectedly reading from line "
            << i << std::endl;
        exit(-1);
      }
      datastream >> m;
      Particle particle(x, v, m);
      particleContainer.addParticle(particle);

      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << std::endl;
    }
  } else {
    std::cout << "Error: could not open file " << filename << std::endl;
    exit(-1);
  }
}

CuboidFileReader::CuboidFileReader() = default;

//CuboidFileReader::~CuboidFileReader() = default;

void CuboidFileReader::readFile(ParticleContainer& particleContainer, char* filename) {
  std::array<double, 3> x;
  std::array<double, 3> v;
  double m;
  std::array<int, 3> s;
  int num_particles = 0;

  double h; //distance
  int e; //Param lennard-jones
  int o; //Param lennard-jones
  double mv; //mean velocity

  std::ifstream input_file(filename);
  std::string tmp_string;

  if (input_file.is_open()) {
    // Read number of particles
    getline(input_file, tmp_string);
    std::istringstream numstream(tmp_string);
    numstream >> num_particles;

    // Skip empty lines and comment lines
    while (tmp_string.empty() || tmp_string[0] == '#') {
      getline(input_file, tmp_string);
    }

    // Read particles
    for (int i = 0; i < num_particles; i++) {
      // Read position and velocity
      for (auto& xj : x) {
        input_file >> xj;
      }
      for (auto& vj : v) {
        input_file >> vj;
      }
      input_file >> m;
      for (auto& vs : s) {
        input_file >> vs;
      }


      // Read additional parameters

      input_file >> h >> e >> o >> mv;

      // Create particle and add to container
      Particle particle(x, v, m);
      ParticleGenerator generator(particle, s.at(0), s.at(1), s.at(2), h, mv);
      particleContainer.addCube(generator.getCube());
    }
    input_file.close();
  } else {
    std::cerr << "Error: could not open file " << filename << std::endl;
    exit(EXIT_FAILURE);
  }
}

