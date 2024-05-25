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
      for (auto particle : particles) {
        pg.generateCuboid(particle, in.cuboids()[i].n1(), in.cuboids()[i].n2(),
                          in.cuboids()[i].n3(), in.cuboids()[i].distance(),
                          in.cuboids()[i].meanVelocity(),
                          in.cuboids()[i].dimension());
      }
      particleContainer.addCube(pg.getCube());
    }
  }
}
