
#include "inputReader/FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"

#include "Calculations.h"
#include "inputReader/XMLReader.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <spdlog/spdlog.h>

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration, std::string outputType, std::string baseName,
                   std::string outputPath, LCParticleContainer& particles);


int main(int argc, char *argsv[]) {
  XMLReader xmlReader(argsv[1]);
  auto start = std::chrono::high_resolution_clock::now();

  std::string particleContainerType = xmlReader.getParticleContainerType();
  std::string inputType = xmlReader.getInputType();
  std::array<double, 3> times = xmlReader.getTime();
  std::string outputType = xmlReader.getOutputType();
  std::string baseName = xmlReader.getBaseName();
  //int writeFrequency = xmlReader.getWriteFrequency();
  xml_schema::boolean performanceMeasurement =
      xmlReader.getPerformanceMeasurement();
  std::string logLevel = xmlReader.getLogLevel();

  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc < 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }
  spdlog::default_logger()->set_level(spdlog::level::info);

  if (logLevel == "trace") {
    spdlog::default_logger()->set_level(spdlog::level::trace);
  } else if (logLevel == "debug") {
    spdlog::default_logger()->set_level(spdlog::level::debug);
  } else if (logLevel == "warn") {
    spdlog::default_logger()->set_level(spdlog::level::warn);
  } else if (logLevel == "error") {
    spdlog::default_logger()->set_level(spdlog::level::err);
  } else if (logLevel == "info") {
    spdlog::default_logger()->set_level(spdlog::level::info);
  } else {
    std::cout << "Invalid log level! " << std::endl;
    return 1;
  }

  double start_time = times[0];
  double end_time = times[1];
  double delta_t = times[2];

  double current_time = start_time;

  int iteration = 0;
  LCParticleContainer lcParticles;
  ParticleContainer normParticles;
  if (particleContainerType == "LC") {
    xmlReader.readXML_LC(lcParticles);
    spdlog::info("read it!");
    spdlog::info("there are {} particles in the Simulation", lcParticles.getParticles().size());
    spdlog::info("there are {} Cells in the Simulation", lcParticles.getCells().size());
  } else {
    xmlReader.readXML(normParticles);
  }

  Calculations normCalculations(normParticles);
  Calculations lcCaluclations(lcParticles);
  spdlog::info("Simulation started with parameters: start_time: {}, end_time: "
               "{}, delta_t: {}, inputType: {}, outputType: {}, baseName: {}, "
               "logLevel: {}, performanceMeasurement: {}, {} "
               "particles, {} cuboids, {} disks ",
               start_time, end_time, delta_t, inputType, outputType, baseName,
               logLevel, performanceMeasurement, lcParticles.getParticles().size(),
               xmlReader.getNumberOfCuboids(), xmlReader.getNumberOfDisks());
  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    if(particleContainerType == "LC") {
      lcCaluclations.calculateX(delta_t);
    } else {
      normCalculations.calculateX(delta_t);
    }

    // calculate new f
    if (inputType == "SF") {
      if(particleContainerType == "LC") {
        lcCaluclations.calculateF();
      } else {
        normCalculations.calculateF();
      }
    } else if(particleContainerType == "LC") {
      lcParticles.handleLJFCalculation();



    } else {
        normCalculations.calculateLJF();
    }

    // calculate new v
    if(particleContainerType == "LC") {
      lcCaluclations.calculateV(delta_t);
    } else {
      normCalculations.calculateV(delta_t);
    }



    iteration++;
    if (!performanceMeasurement) {
      if (iteration % 10 == 0) {
        plotParticles(iteration, outputType, baseName, "../output", lcParticles);
      }
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
void plotParticles(int iteration, std::string outputType, std::string baseName,
                   std::string outputPath, LCParticleContainer& particles) {
  std::filesystem::path dir(outputPath);
  if (!std::filesystem::exists(dir)) {
    std::filesystem::create_directories(dir);
  }
  std::string out_name = outputPath + "/" + baseName;
  if (outputType == "xyz") {

    outputWriter::XYZWriter writer;
    writer.plotParticles(particles, out_name, iteration);
  }

  else {
    outputWriter::VTKWriter vtkWriter;
    int numParticles = particles.getParticles().size();
    vtkWriter.initializeOutput(numParticles);

    for (auto &particle : particles.getParticles()) {
      vtkWriter.plotParticle(particle);
    }

    std::string filename = out_name;
    vtkWriter.writeFile(filename, iteration);
  }
}
