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

#include "../ParticleGenerator.h"

//TODO: sigma und epsilon hinzufügen
FileReader::FileReader() = default;

FileReader::~FileReader() = default;

void FileReader::readFile(ParticleContainer &particleContainer,
                          char *filename) {

  std::array<double, 3> x; //xyz-coord
  std::array<double, 3> v; //velocity
  std::array<double, 3> f; //force
  double m; //mass
  int t; //type
  double epsilon;
  double sigma;
  int num_particles = 0;

  std::ifstream input_file(filename);
  std::string tmp_string;

  if (input_file.is_open()) {

    getline(input_file, tmp_string);
    spdlog::debug("Read line: {}", tmp_string);

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      spdlog::debug("Read line: {}", tmp_string);
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    spdlog::debug("Number of particles: {}", num_particles);

    // Read the next line containing the first particle data
    getline(input_file, tmp_string);
    //spdlog::info("Read line: {}", tmp_string);

    for (int i = 0; i < num_particles; i++) {
      std::istringstream datastream(tmp_string);

      for (auto &xj : x) {
        datastream >> xj;
      }
      for (auto &vj : v) {
        datastream >> vj;
      }

      for (auto &fj : f) {
        datastream >> fj;
      }

      if (datastream.eof()) {
        spdlog::error(
            "Error reading file: eof reached unexpectedly reading from line {}",
            i);
        exit(-1);
      }
      datastream >> m;
      datastream >> t;
      datastream >> epsilon;
      datastream >> sigma;


      Particle particle(x, v, m, t, epsilon, sigma);
      particle.setF(f);
      particleContainer.addParticle(particle);

      getline(input_file, tmp_string);
    }
  } else {
    spdlog::error("Error opening file: {}", filename);
    exit(-1);
  }
}

CuboidFileReader::CuboidFileReader() = default;

// CuboidFileReader::~CuboidFileReader() = default;

void CuboidFileReader::readFileCuboid(ParticleContainer &particleContainer,
                                      char *filename) {
  std::array<double, 3> x;
  std::array<double, 3> v;

  double m;
  int num_particles = 0;
  std::array<double, 3> s; // size im n1 n2 n3 format
  double mv;
  double distance;
  int dimension;
  double temp_init;

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

    // Read the next line containing the first particle data
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
      datastream >> distance;
      datastream >> mv;
      datastream >> dimension;
      if (datastream.eof()) {
        spdlog::error(
            "Error reading file: eof reached unexpectedly reading from line {}",
            i);
        exit(-1);
      }
      datastream >> m;
      datastream >> temp_init;
      Particle particle(x, v, m);
      ParticleGenerator generator;
      if(dimension == 2) {
        generator.generateCuboid(particle, s.at(0), s.at(1), 1, distance,
                               mv, dimension, temp_init);
      } else {
        generator.generateCuboid(particle, s.at(0), s.at(1), s.at(2), distance,
                                  mv, dimension, 0.0);
      }
      particleContainer.addMultipleParticles(generator.getAllParticles());

      getline(input_file, tmp_string);
      if(i != num_particles - 1){
        spdlog::info("Read line: {}", tmp_string);
      }
    }
  } else {
    spdlog::error("Error opening file: {}", filename);
    exit(-1);
  }
}

DiskFileReader::DiskFileReader() = default;

// DiskFileReader::~DiskFileReader() = default;

void DiskFileReader::readFileDisk(ParticleContainer &particleContainer,
                                      char *filename) {
  std::array<double, 3> x;
  std::array<double, 3> v;
  double m;
  int num_particles = 0;
  int radius;
  double distance;
  double mv;
  int dimension;
  double temp_init;

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

    // Read the next line containing the first particle data
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
      datastream >> radius;
      datastream >> distance;
      datastream >> mv;
      datastream >> dimension;
      if (datastream.eof()) {
        spdlog::error(
            "Error reading file: eof reached unexpectedly reading from line {}",
            i);
        exit(-1);
      }
      datastream >> m;
      datastream >> temp_init;
      Particle particle(x, v, m);
      ParticleGenerator generator;
      generator.generateDisk(particle, radius, distance, mv, dimension, temp_init);
      particleContainer.addMultipleParticles(generator.getAllParticles());

      getline(input_file, tmp_string);
      if(i != num_particles - 1){
        spdlog::info("Read line: {}", tmp_string);
      }
    }
  } else {
    spdlog::error("Error opening file: {}", filename);
    exit(-1);
  }
}
