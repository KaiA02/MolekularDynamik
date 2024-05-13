
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h" //Task 3
#include "utils/ArrayUtils.h"

#include <iostream>
#include <list>

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX(double delta_t);

/**
 * calculate the position for all particles
 */
void calculateV(double delta_t);

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration, int fileType);


ParticleContainer particles;

int main(int argc, char *argsv[]) {

  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }



  double start_time = std::stod(argsv[2]);
  double end_time = std::stod(argsv[3]);
  double delta_t = std::stod(argsv[4]);
  int fileType = std::stoi(argsv[5]);
  int fileInputType = std::stoi(argsv[6]);

  if(fileInputType == 1) {
    FileReader fileReader;
    fileReader.readFile(particles, argsv[1]);
  } else {
    CuboidFileReader fileReader;
    fileReader.readFile(particles, argsv[1]);
  }

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calculateX(delta_t);
    // calculate new f
    calculateF();
    // calculate new v
    calculateV(delta_t);

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration, fileType);
    }
    std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}


/**
 * @brief the following function calculates the new forces
 *
 * therefore we use an interator and also for-loops
 * to calculate the physics behind the forces
 */
void calculateF() {

  for (auto &p1 : particles) {
    std::array<double, 3> newForce = {0.0, 0.0, 0.0};

    for (auto &p2 : particles) {
      if (&p1 != &p2) {
        double distSquared = 0.0;
        for (int i = 0; i < 3; ++i) {
          distSquared += std::pow(p1.getX()[i] - p2.getX()[i], 2);
        }

        double distCubed = std::pow(distSquared, 1.5);

        double prefactor = (p1.getM() * p2.getM()) / distCubed;

        for (int i = 0; i < 3; ++i) {
          newForce[i] += prefactor * (p2.getX()[i] - p1.getX()[i]);
        }
      }
    }
    p1.setF(newForce);
  }
}

/**
 * @brief the following function calculates the new positions
 *
 * therefore we use a for-loops
 * to calculate the physics behind the positions (Velocity-Störmer-Verlet)
 */
void calculateX(double delta_t) {
  for (auto &p : particles) {
    std::array<double, 3> newPosition;
    for (int i = 0; i < 3; ++i) {
        newPosition[i] = p.getX()[i] + delta_t * p.getV()[i] + delta_t * delta_t * (p.getF()[i] / (2 * p.getM()));
    }
    p.setX(newPosition);
  }
}

/**
 * @brief the following function calculates the new velocities
 *
 * therefore we use a for-loops
 * to calculate the physics behind the velocities (Velocity-Störmer-Verlet)
 */
void calculateV(double delta_t) {
  for (auto &p : particles) {
    std::array<double, 3> newVelocity;
    for (int i = 0; i < 3; ++i) {
      newVelocity[i] = p.getV()[i] + delta_t * (p.getF()[i] + p.getOldF()[i]) / (2 * p.getM());
    }
    p.setV(newVelocity);
  }
}

/**
 * @brief the following function plots the particles
 *
 * in this funnction the particles are plotted in xyz files and also in vtu files
 *
 * @param iteration is the number of iterations of the particles
 */
void plotParticles(int iteration, int fileType) {

  if(fileType == 2) {
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
