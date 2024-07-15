
#include "Calculations.h"
#include "Container/LCParticleContainer.h"
#include "Particle.h"
#include "inputReader/FileReader.h"
#include "inputReader/XMLReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <spdlog/spdlog.h>

#include "Thermostat.h"

/**
 * plot the particles to a xyz-file
 */
void plotParticlesLC(int iteration, std::string outputType,
                     std::string baseName, std::string outputPath,
                     LCParticleContainer &particles);
void plotParticles(int iteration, std::string outputType, std::string baseName,
                   std::string outputPath, ParticleContainer &particles);
void displayProgressBar(int progress, int total,
                        std::chrono::high_resolution_clock::time_point start);

void saveState(std::vector<Particle> particles);

void outputDiffusionToFile(const std::vector<std::pair<double, int>> &diffusion,
                           const std::string &filename, std::string outputPath);
void rdfOutput(std::vector<std::pair<int, std::map<double, double>>> rdf,
               std::string filename, std::string outputPath);

int main(int argc, char *argsv[]) {

  LCParticleContainer lcParticles;
  ParticleContainer normParticles;

  XMLReader xmlReader(argsv[1]);

  std::string particleContainerType = xmlReader.getParticleContainerType();
  bool smoothLJ = xmlReader.getSmoothLJ();
  double r_l = xmlReader.getR_L();
  std::array<double, 3> times = xmlReader.getTime();
  std::string outputType = xmlReader.getOutputType();
  std::string baseName = xmlReader.getBaseName();
  // int writeFrequency = xmlReader.getWriteFrequency();
  xml_schema::boolean performanceMeasurement =
      xmlReader.getPerformanceMeasurement();
  std::string logLevel = xmlReader.getLogLevel();

  std::string thermostatOn = xmlReader.ThermostatON();
  double temp_init = xmlReader.getTemp_init();
  int n_thermostat = xmlReader.getN_Thermostat();
  double temp_target = xmlReader.getTemp_Target();
  double delta_temp = xmlReader.getDelta_Temp();
  bool statisticsOn = xmlReader.getStatisticsOn();
  double rdfDeltaR = xmlReader.getRdfDeltaR();
  Thermostat thermostat(temp_init, n_thermostat, temp_target, delta_temp);
  std::vector<Particle> prevParticles;
  std::vector<std::pair<double, int>> diffusion;
  std::vector<std::pair<int, std::map<double, double>>> rdf;

  // std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc > 2) {
    FileReader stateReader;
    stateReader.readFile(lcParticles, argsv[2]);
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

  double totalIterations = (end_time - start_time) / delta_t;
  spdlog::debug("MolSim: end_time: {}", end_time);
  spdlog::debug("MolSim: start_time: {}", start_time);
  spdlog::debug("MolSim: delta_t: {}", delta_t);
  int progress = 0;

  int molecule_updates = 0;

  int iteration = 0;

  if (particleContainerType == "LC") {
    xmlReader.readXML_LC(lcParticles);
  } else {
    xmlReader.readXML(normParticles);
  }
  spdlog::info("read");
  Calculations normCalculations(normParticles);
  Calculations lcCaluclations(lcParticles);
  lcCaluclations.setR_cutoff(lcParticles.getR_cutoff());
  lcCaluclations.setG_grav(lcParticles.getG_grav());
  lcCaluclations.setSmoothLJ(smoothLJ);
  lcCaluclations.setR_L(r_l);
  lcParticles.setUpEpsilonAndSigmas();
  spdlog::warn("Simulation started with parameters: start_time: {}, end_time: "
               "{}, delta_t: {}, temp_init: {}, temp_target: {}, smoothLJ: "
               "{}, outputType: {}, baseName: {}, "
               "logLevel: {}, performanceMeasurement: {}, {} "
               "particles, {} cuboids, {} disks ",
               start_time, end_time, delta_t, temp_init, temp_target, smoothLJ,
               outputType, baseName, logLevel, performanceMeasurement,
               lcParticles.getParticles().size(),
               xmlReader.getNumberOfCuboids(), xmlReader.getNumberOfDisks());
  auto start = std::chrono::high_resolution_clock::now();
  // for this loop, we assume: current x, current f and current v are known
  if (particleContainerType == "LC") {

    if (thermostatOn == "YES") {
      // temperature setting
      spdlog::info("current Temperature: {}",
                   thermostat.getCurrentTemp(lcParticles.getParticles()));
      thermostat.setInitialTemperature(lcParticles.getParticles());
      spdlog::info("current Temperature: {}",
                   thermostat.getCurrentTemp(lcParticles.getParticles()));
      molecule_updates += lcParticles.getParticles().size();
    }
    prevParticles = lcParticles.getParticles();

    while (current_time < end_time) {

      lcCaluclations.calculateX(delta_t);

      molecule_updates += lcParticles.getParticles().size();
      lcParticles.handleLJFCalculation(lcCaluclations, int(current_time));

      molecule_updates +=
          5 * lcParticles.getParticles().size(); // only provisionally
      lcCaluclations.calculateV(delta_t);

      molecule_updates += lcParticles.getParticles().size();
      iteration++;

      if (thermostatOn == "YES") {
        if (n_thermostat == 0) {
          thermostat.gradualScaling(lcParticles.getParticles());
        } else {
          if (iteration % n_thermostat == 0) {
            if (delta_temp == 0) {
              thermostat.setTemperatureDirectly(lcParticles.getParticles());
            } else {
              thermostat.gradualScaling(lcParticles.getParticles());
            }
          }
        }
        molecule_updates += lcParticles.getParticles().size();
      }
      spdlog::debug("MolSim: applied Thermostat");
      if (!performanceMeasurement) {
        if (iteration % 10 == 0) {
          if (iteration % 1000 == 0) {
            diffusion.emplace_back(
                Calculations::calculateDiffusion(lcParticles.getParticles(),
                                                 prevParticles),
                iteration);
            rdf.emplace_back(iteration,
                             lcCaluclations.calculateLocalDensities(
                                 lcParticles.getParticles(), rdfDeltaR));
            for (const auto &particle : lcParticles.getParticles()) {
              prevParticles.push_back(particle);
            }
          }
          if (thermostatOn == "YES") {
            spdlog::debug(
                "current Temperature: {}",
                thermostat.getCurrentTemp(lcParticles.getParticles()));
          }
          plotParticlesLC(iteration, outputType, baseName, "../output",
                          lcParticles);
          displayProgressBar(progress, totalIterations, start);
        }
      }
      progress++;
      current_time += delta_t;
    }
    saveState(lcParticles.getParticles());
  } else {
    while (current_time < end_time) {
      normCalculations.calculateX(delta_t);
      normCalculations.calculateF();
      normCalculations.calculateV(delta_t);
      iteration++;

      if (!performanceMeasurement) {
        if (iteration % 10 == 0) {
          if (iteration % 1000 == 0) {
            diffusion.emplace_back(
                Calculations::calculateDiffusion(lcParticles.getParticles(),
                                                 prevParticles),
                iteration);
            rdf.emplace_back(iteration,
                             lcCaluclations.calculateLocalDensities(
                                 lcParticles.getParticles(), rdfDeltaR));
            prevParticles = lcParticles.getParticles();
          }
          plotParticles(iteration, outputType, baseName, "../output",
                        normParticles);
        }
      }
      current_time += delta_t;
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  molecule_updates = molecule_updates / elapsed.count();
  spdlog::warn("Simulation finished in {} seconds", elapsed.count());
  if (performanceMeasurement) {
    spdlog::warn("There were {} molecule-updates per second", molecule_updates);
  }
  int parkedCounter = 0;
  for (auto p : lcParticles.getParticles()) {
    if (p.getX().at(0) < -5.0) {
      parkedCounter++;
    }
  }
  spdlog::warn("there are {} parked Particles", parkedCounter);
  /**for (int m=0;m<diffusion.size();m++){
    std::cout << "Diffusion at iteration " << diffusion[m].second << " is " <<
  diffusion[m].first << std::endl;
  }**/
  outputDiffusionToFile(diffusion, "diffusion.txt", "../output/statistics");
  rdfOutput(rdf, "rdf.txt", "../output/statistics");
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
void plotParticlesLC(int iteration, std::string outputType,
                     std::string baseName, std::string outputPath,
                     LCParticleContainer &particles) {
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
void plotParticles(int iteration, std::string outputType, std::string baseName,
                   std::string outputPath, ParticleContainer &particles) {
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

void displayProgressBar(int progress, int total,
                        std::chrono::high_resolution_clock::time_point start) {
  const int barWidth = 100;

  std::cout << "[";
  int pos = barWidth * progress / total;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)
      std::cout << "=";
    else if (i == pos)
      std::cout << ">";
    else
      std::cout << " ";
  }

  auto now = std::chrono::high_resolution_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
  auto totalEstimatedTime = total > 0 ? (elapsed * total) / progress : 0;
  auto remainingTime = totalEstimatedTime - elapsed;

  std::cout << "] " << int(progress * 100.0 / total)
            << " %, estimated time remaining: " << remainingTime << "s\r";
  std::cout.flush();
}

void saveState(std::vector<Particle> particles) {
  std::ofstream outFile("../input/output.txt");

  // Write the number of particles to the file
  outFile << particles.size() << "\n";

  // Iterate over the particles and write their attributes to the file
  for (const auto &particle : particles) {
    // Write the xyz-coordinates
    for (const auto &coord : particle.getX()) {
      outFile << coord << " ";
    }

    // Write the velocities
    for (const auto &velocity : particle.getV()) {
      outFile << velocity << " ";
    }

    // Write the force
    for (const auto &force : particle.getF()) {
      outFile << force << " ";
    }

    // Write the mass
    outFile << particle.getM() << " ";

    // Write the type
    outFile << particle.getType() << " ";

    // Write the epsilon
    outFile << particle.getEpsilon() << " ";

    // Write the sigma
    outFile << particle.getSigma() << "\n";
  }

  outFile.close();
  spdlog::warn("State is saved");
}

void rdfOutput(std::vector<std::pair<int, std::map<double, double>>> rdf,
               std::string filename, std::string outputPath) {
  std::filesystem::path dir(outputPath);
  if (!std::filesystem::exists(dir)) {
    std::filesystem::create_directories(dir);
  }
  std::string out_name = outputPath + "/" + filename;
  std::ofstream file(out_name);
  if (!file.is_open()) {
    std::cerr << "Failed to open file for writing: " << out_name << std::endl;
    return;
  }

  for (const auto &pair : rdf) {
    file << pair.first << " ";
    for (const auto &innerPair : pair.second) {
      file << innerPair.first << " " << innerPair.second << " ";
    }
    file << "\n";
  }
  file.close();
}

void outputDiffusionToFile(const std::vector<std::pair<double, int>> &diffusion,
                           const std::string &filename,
                           std::string outputPath) {
  std::filesystem::path dir(outputPath);
  if (!std::filesystem::exists(dir)) {
    std::filesystem::create_directories(dir);
  }
  std::string out_name = outputPath + "/" + filename;
  std::ofstream file(out_name);
  if (!file.is_open()) {
    std::cerr << "Failed to open file for writing: " << out_name << std::endl;
    return;
  }

  for (const auto &pair : diffusion) {
    file << pair.second << " " << pair.first
         << "\n"; // iteration number and diffusion value
  }

  file.close();
}