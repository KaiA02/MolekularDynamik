//
// Created by jh on 18.06.2024.
//

#include "OutputToInput.h"

#include "spdlog/spdlog.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>


OutputToInput::OutputToInput() = default;

OutputToInput::~OutputToInput() = default;

void OutputToInput::readOutput(LCParticleContainer &particleContainer,
                          char *filename) {

    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    int num_particles = 0;
    int type = 0;

    std::ifstream input_file(filename);
    std::string tmp_string;

    if (input_file.is_open()) {

        getline(input_file, tmp_string);
        spdlog::info("Read line: {}", tmp_string);

        while (tmp_string.empty() or tmp_string[0] == '#') {
            getline(input_file, tmp_string);
            spdlog::info("Read line: {}", tmp_string);
        }

        if(startsWith(tmp_string, "<Piece")) {
            getline(input_file, tmp_string);
            num_particles = extractNumberOfPoints(tmp_string);
        }


        spdlog::info("Number of particles: {}", num_particles);

        // Read the next line containing the first particle data
        getline(input_file, tmp_string);
        spdlog::info("Read line: {}", tmp_string);

        for (int i = 0; i < num_particles; i++) {
            std::istringstream datastream(tmp_string);

            for (auto &xj : x) {
                datastream >> xj;
            }
            for (auto &vj : v) {
                datastream >> vj;
            }

            if (datastream.eof()) {
                spdlog::error(
                    "Error reading file: eof reached unexpectedly reading from line {}",
                    i);
                exit(-1);
            }
            datastream >> m;

            Particle particle(x, v, m);
            particleContainer.addParticle(particle);

            getline(input_file, tmp_string);
            if(i != num_particles - 1){
                spdlog::info("Read line: {}", tmp_string);
            }
        }
    } else {
        spdlog::error("Error opening file: {}", filename);
        exit(-1);
    }
}

bool OutputToInput::startsWith(const std::string &str, const std::string &prefix) {
    return str.compare(0, prefix.size(), prefix) == 0;
}

int OutputToInput::extractNumberOfPoints(const std::string &line) {
    std::string keyword = "NumberOfPoints=\"";
    size_t startPos = line.find(keyword);
    if (startPos != std::string::npos) {
        startPos += keyword.length();
        size_t endPos = line.find("\"", startPos);
        if (endPos != std::string::npos) {
            std::string numberStr = line.substr(startPos, endPos - startPos);
            return std::stoi(numberStr);
        }
    }
    throw std::runtime_error("NumberOfPoints not found in line");
}


