
#include "inputReader/FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "Container/LCParticleContainer.h"
#include "Calculations.h"
#include "inputReader/XMLReader.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <spdlog/spdlog.h>

#include "Thermostat.h"

/**
 * plot the particles to a xyz-file
 */
void plotParticlesLC(int iteration, std::string outputType, std::string baseName,
                   std::string outputPath, LCParticleContainer& particles);
void plotParticles(int iteration, std::string outputType, std::string baseName,
                   std::string outputPath, ParticleContainer& particles);
void displayProgressBar(int progress, int total, std::chrono::high_resolution_clock::time_point start);


int main(int argc, char *argsv[]) {
  XMLReader xmlReader(argsv[1]);

  std::string particleContainerType = xmlReader.getParticleContainerType();
  std::string inputType = xmlReader.getInputType();
  std::array<double, 3> times = xmlReader.getTime();
  std::string outputType = xmlReader.getOutputType();
  std::string baseName = xmlReader.getBaseName();
  //int writeFrequency = xmlReader.getWriteFrequency();
  xml_schema::boolean performanceMeasurement =
      xmlReader.getPerformanceMeasurement();
  std::string logLevel = xmlReader.getLogLevel();

  std::string thermostatOn = xmlReader.ThermostatON();
  double temp_init = xmlReader.getTemp_init();
  int n_thermostat = xmlReader.getN_Thermostat();
  double temp_target = xmlReader.getTemp_Target();
  double delta_temp = xmlReader.getDelta_Temp();
  Thermostat thermostat(temp_init, n_thermostat, temp_target, delta_temp);

  //std::cout << "Hello from MolSim for PSE!" << std::endl;
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

  double totalIterations = (end_time - start_time) / delta_t;
  int progress = 0;

  int iteration = 0;
  LCParticleContainer lcParticles;
  ParticleContainer normParticles;
  if (particleContainerType == "LC") {
    xmlReader.readXML_LC(lcParticles);
  } else {
    xmlReader.readXML(normParticles);
  }

  Calculations normCalculations(normParticles);
  Calculations lcCaluclations(lcParticles);
  lcCaluclations.setR_cutoff(lcParticles.getR_cutoff());
  lcCaluclations.setG_grav(lcParticles.getG_grav());
  spdlog::warn("Simulation started with parameters: start_time: {}, end_time: "
               "{}, delta_t: {}, temp_init: {}, temp_target: {}, inputType: {}, outputType: {}, baseName: {}, "
               "logLevel: {}, performanceMeasurement: {}, {} "
               "particles, {} cuboids, {} disks ",
               start_time, end_time, delta_t, temp_init, temp_target, inputType, outputType, baseName,
               logLevel, performanceMeasurement, lcParticles.getParticles().size(),
               xmlReader.getNumberOfCuboids(), xmlReader.getNumberOfDisks());
  auto start = std::chrono::high_resolution_clock::now();
  // for this loop, we assume: current x, current f and current v are known
  if(particleContainerType == "LC") {

    if(thermostatOn == "YES") {
    //temperature setting
     spdlog::info("current Temperature: {}", thermostat.getCurrentTemp(lcParticles.getParticles()));
     thermostat.setInitialTemperature(lcParticles.getParticles());
     spdlog::info("current Temperature: {}", thermostat.getCurrentTemp(lcParticles.getParticles()));
    }

    while(current_time < end_time) {
      lcCaluclations.calculateX(delta_t);
      if(inputType == "SF") { //Simple Force
        lcCaluclations.calculateLJF();
      } else {
        lcParticles.handleLJFCalculation(lcCaluclations);
      }


      lcCaluclations.calculateV(delta_t);
      iteration++;

      if(thermostatOn == "YES") {
        if(n_thermostat == 0) {
          thermostat.gradualScaling(lcParticles.getParticles());
        } else {
          if(iteration % n_thermostat == 0) {
            if(delta_temp == 0) {
              thermostat.setTemperatureDirectly(lcParticles.getParticles());
           } else {
              thermostat.gradualScaling(lcParticles.getParticles());
           }
          }
        }
      }

      if (!performanceMeasurement) {
        if (iteration % 10 == 0) {
          if(thermostatOn == "YES") {
            spdlog::info("current Temperature: {}", thermostat.getCurrentTemp(lcParticles.getParticles()));
          }
          plotParticlesLC(iteration, outputType, baseName, "../output", lcParticles);
          displayProgressBar(progress, totalIterations, start);
        }
      }
      progress++;
      current_time += delta_t;

    }
  } else {

    if(thermostatOn == "YES") {
      //temperature setting
      spdlog::info("current Temperature: {}", thermostat.getCurrentTemp(normParticles.getParticles()));
      thermostat.setInitialTemperature(normParticles.getParticles());
      spdlog::info("current Temperature: {}", thermostat.getCurrentTemp(normParticles.getParticles()));
    }

    while(current_time < end_time) {
      normCalculations.calculateX(delta_t);
      if(inputType == "SF") { //Simple Force
        normCalculations.calculateF();
      } else {
        normCalculations.calculateLJF();
      }
      normCalculations.calculateV(delta_t);
      iteration++;

      if(thermostatOn == "YES") {
        if(n_thermostat == 0) {
         thermostat.gradualScaling(normParticles.getParticles());
        } else {
          if(iteration % n_thermostat == 0) {
           if(delta_temp == 0) {
             thermostat.setTemperatureDirectly(normParticles.getParticles());
           } else {
             thermostat.gradualScaling(normParticles.getParticles());
           }
         }
        }
      }

      if (!performanceMeasurement) {
        if (iteration % 10 == 0) {
          if(thermostatOn == "YES") {
            spdlog::info("current Temperature: {}", thermostat.getCurrentTemp(normParticles.getParticles()));
          }
          plotParticles(iteration, outputType, baseName, "../output", normParticles);
        }
      }
      current_time += delta_t;
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  spdlog::warn("Simulation finished in {} seconds", elapsed.count());
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
void plotParticlesLC(int iteration, std::string outputType, std::string baseName,
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
void plotParticles(int iteration, std::string outputType, std::string baseName,
                   std::string outputPath, ParticleContainer& particles) {
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

void displayProgressBar(int progress, int total, std::chrono::high_resolution_clock::time_point start) {
  const int barWidth = 100;

  std::cout << "[";
  int pos = barWidth * progress / total;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }

  auto now = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
  auto totalEstimatedTime = total > 0 ? (elapsed * total) / progress : 0;
  auto remainingTime = totalEstimatedTime - elapsed;

  std::cout << "] " << int(progress * 100.0 / total) << " %, estimated time remaining: " << remainingTime << "s\r";
  std::cout.flush();
}
