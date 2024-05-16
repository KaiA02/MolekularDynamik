/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"

#include "spdlog/spdlog.h"

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
    spdlog::info("Read line: {}", tmp_string);

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      spdlog::info("Read line: {}", tmp_string);
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    spdlog::info("Number of particles: {}", num_particles);
    getline(input_file, tmp_string);
    spdlog::info("Read line: {}", tmp_string);

    for (int i = 0; i < num_particles; i++) {
      std::istringstream datastream(tmp_string);

      for (auto &xj : x) {
        datastream >> xj;
      }
      for (auto &vj : v) {
        datastream >> vj;
      }
      if (datastream.eof()) {
        spdlog::error(
            "Error reading file: eof reached unexpectedly reading from line {}",
            i);
        exit(-1);
      }
      datastream >> m;
      Particle particle(x, v, m);
      particleContainer.addParticle(particle);

      getline(input_file, tmp_string);
      spdlog::info("Read line: {}", tmp_string);
    }
  } else {
    spdlog::error("Error opening file: {}", filename);
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
      spdlog::info("Read line: {}", tmp_string);

      while (tmp_string.empty() or tmp_string[0] == '#') {
        getline(input_file, tmp_string);
        spdlog::info("Read line: {}", tmp_string);
      }

      std::istringstream numstream(tmp_string);
      numstream >> num_particles;
      //numstream >> h;
      //numstream >> mv;

      spdlog::info("Number of particles: {}", num_particles);
      getline(input_file, tmp_string);
      spdlog::info("Read line: {}", tmp_string);

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
          spdlog::error(
             "Error reading file: eof reached unexpectedly reading from line {}",
             i);
          exit(-1);
        }
        datastream >> m;
        Particle particle(x, v, m);
        ParticleGenerator generator;
        generator.generateCuboid(particle, s.at(0), s.at(1), s.at(2), 1.1225, 0.1);
        particleContainer.addCube(generator.getCube());

        getline(input_file, tmp_string);
        spdlog::info("Read line: {}", tmp_string);
      }
    } else {
      spdlog::error("Error opening file: {}", filename);
      exit(-1);
    }
  }

