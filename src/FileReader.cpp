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

void CuboidFileReader::readFileCuboid(ParticleContainer& particleContainer, char* filename) {
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    int num_particles = 0;
    //double h; //distance
    std::array<double, 3> s; //size im n1 n2 n3 format
    //double mv; //meanVelocityInput


    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {

      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << std::endl;

      while (tmp_string.empty() or tmp_string[0] == '#') {
        getline(input_file, tmp_string);
        //std::cout << "Read line: " << tmp_string << std::endl;
      }

      std::istringstream numstream(tmp_string);
      numstream >> num_particles;
      //numstream >> h;
      //numstream >> mv;

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
        for (auto &sj : s) {
          datastream >> sj;
        }
        if (datastream.eof()) {
          std::cout
              << "Error reading file: eof reached unexpectedly reading from line "
              << i << std::endl;
          exit(-1);
        }
        datastream >> m;
        Particle particle(x, v, m);
        std::cout << "Particle: " << particle.toString() << std::endl;
        ParticleGenerator generator(particle, s.at(0), s.at(1), s.at(2), 1.1225, 0.1);
        particleContainer.addCube(generator.getCube());

        getline(input_file, tmp_string);
        std::cout << "Read line: " << tmp_string << std::endl;
      }
    } else {
      std::cout << "Error: could not open file " << filename << std::endl;
      exit(-1);
    }
  }

