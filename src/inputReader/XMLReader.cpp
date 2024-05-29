//
// Created by joshu on 25.05.2024.
//
#include "XMLReader.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "simulation.hxx"
#include <fstream>
#include <iostream>
#include <memory>

XMLReader::XMLReader(const std::string &filePath) {
  try {
    sim = simulation_(filePath);
  } catch (const xml_schema::exception &e) {
    std::cerr << e << std::endl;
    throw;
  }
}
void XMLReader::readXML(ParticleContainer &particleContainer) {
  input in = sim->input();
  std::vector<Particle> particles;
  for (int i = 0; i < in.particles().size(); i++) {
    Particle p(
        {in.particles()[i].x(), in.particles()[i].y(), in.particles()[i].z()},
        {in.particles()[i].velocityX(), in.particles()[i].velocityY(),
         in.particles()[i].velocityZ()},
        in.particles()[i].mass());
    particleContainer.addParticle(p);
    particles.push_back(p);
  }
  if (in.cuboids().size() > 0) {
    std::vector<std::vector<Particle>> cubes;
    particleContainer.resetParticles();
    for (int i = 0; i < in.cuboids().size(); i++) {
      ParticleGenerator pg;
      Particle particle = particles[i];
      int dimension = in.cuboids()[i].dimension();
      int n3 = in.cuboids()[i].n3();
      if (dimension == 2) {
        n3 = 1;
      }
      pg.generateCuboid(particle, in.cuboids()[i].n1(), in.cuboids()[i].n2(),
                        n3, in.cuboids()[i].distance(),
                        in.cuboids()[i].meanVelocity(), dimension);

      particleContainer.addCube(pg.getCube());
    }
  }
}

std::array<double, 3> XMLReader::getTime() {
  std::array<double, 3> time;
  time[0] = sim->input().tStart();
  time[1] = sim->input().tEnd();
  time[2] = sim->input().deltaT();
  return time;
}

std::string XMLReader::getInputType() { return sim->input().inputType(); }

std::string XMLReader::getOutputType() { return sim->output().outputType(); }

std::string XMLReader::getBaseName() { return sim->output().baseName(); }

int XMLReader::getWriteFrequency() { return sim->output().writeFrequency(); }

xml_schema::boolean XMLReader::getPerformanceMeasurement() {
  if (sim->config().performanceMeasurement().present()) {
    return sim->config().performanceMeasurement().get();
  }
  return false;
}

std::string XMLReader::getLogLevel() {
  if (sim->config().logLevel().present()) {
    return sim->config().logLevel().get();
  }
  return "info";
}

int XMLReader::getNumberOfParticles() {
  return sim->input().particles().size();
}
int XMLReader::getNumberOfCuboids() { return sim->input().cuboids().size(); }
// int XMLReader::getNumberOfSpheres(){ return sim->input().spheres().size();}
