
#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"

#include "Calculations.h"
#include <chrono>
#include <iostream>
#include <spdlog/spdlog.h>

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration, int fileType);

ParticleContainer particles;
Calculations calculations(particles);

int main(int argc, char *argsv[]) {
  auto start = std::chrono::high_resolution_clock::now();
  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc < 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }
  spdlog::default_logger()->set_level(spdlog::level::info);
  if (argc > 8) {
    std::string log_level = argsv[8];
    if (log_level == "trace") {
      spdlog::default_logger()->set_level(spdlog::level::trace);
    } else if (log_level == "debug") {
      spdlog::default_logger()->set_level(spdlog::level::debug);
    } else if (log_level == "warn") {
      spdlog::default_logger()->set_level(spdlog::level::warn);
    } else if (log_level == "error") {
      spdlog::default_logger()->set_level(spdlog::level::err);
    } else if (log_level == "info") {
      spdlog::default_logger()->set_level(spdlog::level::info);
    } else {
      std::cout << "Invalid log level! " << std::endl;
      return 1;
    }
  }

  double start_time = std::stod(argsv[2]);
  double end_time = std::stod(argsv[3]);
  double delta_t = std::stod(argsv[4]);
  int fileType = std::stoi(argsv[5]);
  int fileInputType = std::stoi(argsv[6]);
  int performanceMeasurement = std::stoi(argsv[7]);

  if (fileInputType == 1) {
    FileReader fileReader;
    fileReader.readFile(particles, argsv[1]);
  } else {
    CuboidFileReader fileReader;
    fileReader.readFileCuboid(particles, argsv[1]);
  }

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calculations.calculateX(delta_t);
    // calculate new f
    if (fileInputType == 1) {
      calculations.calculateF();
    } else {
      calculations.calculateLJF();
    }
    // calculate new v
    calculations.calculateV(delta_t);

    iteration++;
    if (performanceMeasurement != 1) {
      if (iteration % 10 == 0) {
        plotParticles(iteration, fileType);
      }
      spdlog::trace("Iteration {} finished", iteration);
    }

    current_time += delta_t;
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  spdlog::info("Simulation finished in {} seconds", elapsed.count());
  return 0;
}

/**
 * @brief the following function plots the particles
 *
 * in this funnction the particles are plotted in xyz files and also in vtu
 * files
 *
 * @param iteration is the number of iterations of the particles
 */
void plotParticles(int iteration, int fileType) {

  if (fileType == 2) {
    std::string out_name("output_xyz");

    outputWriter::XYZWriter writer;
    writer.plotParticles(particles, out_name, iteration);
  }

  else {
    outputWriter::VTKWriter vtkWriter;
    int numParticles = particles.size();
    vtkWriter.initializeOutput(numParticles);

    for (auto &particle : particles) {
      vtkWriter.plotParticle(particle);
    }

    std::string filename = "MD_vtk";
    vtkWriter.writeFile(filename, iteration);
  }
}
