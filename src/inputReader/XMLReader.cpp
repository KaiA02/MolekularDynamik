//
// Created by joshu on 25.05.2024.
//
#include "XMLReader.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "simulation.hxx"
#include "spdlog/spdlog.h"

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
  for (size_t i = 0; i < in.particles().size(); i++) {
    ParticleGenerator pg;
    Particle particle(
        {in.particles()[i].x(), in.particles()[i].y(), in.particles()[i].z()},
        {in.particles()[i].velocityX(), in.particles()[i].velocityY(),
         in.particles()[i].velocityZ()},
        in.particles()[i].mass());
    particle.setType(i);
    if (i < in.cuboids().size()) {
      int dimension = in.cuboids()[i].dimension();
      int n3 = in.cuboids()[i].n3();
      if (dimension == 2) {
        n3 = 1;
      }
      spdlog::debug("generating Cube");
      pg.generateCuboid(particle, in.cuboids()[i].n1(), in.cuboids()[i].n2(),
                        n3, in.cuboids()[i].distance(),
                        in.cuboids()[i].meanVelocity(), dimension);

      particleContainer.addMultipleParticles(pg.getAllParticles());
    } else if (i < in.disk().size()) {
      spdlog::debug("generating Disk");
      int dimension = in.disk()[i].dimension();
      pg.generateDisk(particle, in.disk()[i].radius(), in.disk()[i].distance(),
                      dimension);
      particleContainer.addMultipleParticles(pg.getAllParticles());
    } else {
      particleContainer.addParticle(particle);

    }
  }

}
void XMLReader::readXML_LC(LCParticleContainer &particleContainer) {
  input in = sim->input();
  particleContainer.generateCells(in.domainSizeX(), in.domainSizeY(), in.domainSizeZ(), in.r_cutoff());
  for (size_t i = 0; i < in.particles().size(); i++) {
    ParticleGenerator pg;
    Particle particle(
        {in.particles()[i].x(), in.particles()[i].y(), in.particles()[i].z()},
        {in.particles()[i].velocityX(), in.particles()[i].velocityY(),
         in.particles()[i].velocityZ()},
        in.particles()[i].mass());
    particle.setType(i);
    if (i < in.cuboids().size()) { //case its a cube
      int dimension = in.cuboids()[i].dimension();
      int n3 = in.cuboids()[i].n3();
      if (dimension == 2) {
        n3 = 1;
      }
      spdlog::debug("generating Cube");
      pg.generateCuboid(particle, in.cuboids()[i].n1(), in.cuboids()[i].n2(),
                        n3, in.cuboids()[i].distance(),
                        in.cuboids()[i].meanVelocity(), dimension);

      particleContainer.addMultipleParticles(pg.getAllParticles());
      spdlog::info("added {} particles to the generator", pg.getAllParticles().size());
      spdlog::info("there are {} particles in the container now", particleContainer.getParticles().size());
    } else if (i < in.disk().size()) { //case its a disk
      //spdlog::debug("generating Disk");
      int dimension = in.disk()[i].dimension();
      pg.generateDisk(particle, in.disk()[i].radius(), in.disk()[i].distance(),
                      dimension);
      particleContainer.addMultipleParticles(pg.getAllParticles());
      //spdlog::info("added {} particles to the container", particleContainer.getParticles().size());
    } else { //case its a single particle
      particleContainer.addParticle(particle);
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

std::string XMLReader::getParticleContainerType() { return sim->input().particleContainerType(); }

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

int XMLReader::getNumberOfDisks() { return sim->input().disk().size(); }
// int XMLReader::getNumberOfSpheres(){ return sim->input().spheres().size();}