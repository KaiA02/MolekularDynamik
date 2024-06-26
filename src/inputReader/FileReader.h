/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "../Particle.h"

#include <list>

#include "../Container/ParticleContainer.h"
/**
 * @brief FileReader class
 */
class FileReader {

public:
  FileReader();
  virtual ~FileReader();

  /**
   * reads a file and adds the resulting particles to the particleContainer
   * @param particleContainer the container to which the particles will be added
   * @param filename the name of the input file
   */
  void readFile(ParticleContainer &particleContainer, char *filename);


};
/**
 * @brief CuboidFileReader extends FileReader
 */
class CuboidFileReader : public FileReader {
public:
  CuboidFileReader();
  // virtual ~CuboidFileReader();
  /**
   * @brief reads a file and adds the resulting cuboid to the
   * particleContainer
   * @param particleContainer
   * @param filename
   */
  void readFileCuboid(ParticleContainer &particleContainer, char *filename);
};

/**
 * @brief DiskFileReader extends FileReader
 */
class DiskFileReader : public FileReader {
public:
 DiskFileReader();
 // virtual ~DiskFileReader();
 /**
  * @brief reads a file and adds the resulting disk to the
  * particleContainer
  * @param particleContainer
  * @param filename
  */
 void readFileDisk(ParticleContainer &particleContainer, char *filename);
};