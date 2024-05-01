
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
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

constexpr double start_time = 0;
constexpr double end_time = 1000;
constexpr double delta_t = 0.014;

// TODO: what data structure to pick?
std::list<Particle> particles;

int main(int argc, char *argsv[]) {

  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
    }
    std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}

void calculateF() {
  std::list<Particle>::iterator iterator;
  iterator = particles.begin();

  for (auto &p1 : particles) {
    std::array<double, 3> newForce = {0.0, 0.0, 0.0};

    for (auto &p2 : particles) {
      // @DONE: insert calculation of forces here!
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

void calculateX() {
  for (auto &p : particles) {
    // @DONE: insert calculation of position updates here!
    std::array<double, 3> newPosition;
    for (int i = 0; i < 3; ++i) {
        newPosition[i] = p.getX()[i] + delta_t * p.getV()[i] + delta_t * delta_t * (p.getF()[i] / (2 * p.getM()));
    }
    p.setX(newPosition);
  }
}

void calculateV() {
  for (auto &p : particles) {
    // @DONE: insert calculation of veclocity updates here!
    std::array<double, 3> newVelocity;
    for (int i = 0; i < 3; ++i) {
      newVelocity[i] = p.getV()[i] + delta_t * (p.getF()[i] + p.getOldF()[i]) / (2 * p.getM()); //getOldF() ist hier vielleicht falsch
    }
    p.setV(newVelocity);
  }
}

void plotParticles(int iteration) {

  std::string out_name("MD_vtk");

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);

  //Task 3
  outputWriter::VTKWriter vtkWriter;
  int numParticles = particles.size();
  vtkWriter.initializeOutput(numParticles);

  for (auto &particle : particles) {
    vtkWriter.plotParticle(particle);
  }

  std::string filename = "output";
  vtkWriter.writeFile(filename, iteration);
}
