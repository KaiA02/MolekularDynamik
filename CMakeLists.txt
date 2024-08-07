cmake_minimum_required(VERSION 3.10)

# define project name, version
project(PSEMolDyn.GroupB VERSION 0.0.1)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Add the -pg flag for gprof profiling
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

# create make target
add_executable(MolSim ${MY_SRC}
        src/ParticleGenerator.cpp
        src/ParticleGenerator.h
        src/Calculations.cpp
        src/Calculations.h
        src/inputReader/XMLReader.cpp
        src/inputReader/XMLReader.h
        src/inputReader/simulation.cxx
        src/inputReader/simulation.hxx
        src/Thermostat.cpp
        src/Thermostat.h
)

# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
        cxx_std_17
)

target_include_directories(MolSim
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

find_package(XercesC 3.2 REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${XercesC_INCLUDE_DIR})

if(OpenMP_CXX_FOUND)
    target_link_libraries(MolSim PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(MolSim
        PUBLIC
        XercesC::XercesC
        spdlog::spdlog
)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(doxygen)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.8.5
)

FetchContent_MakeAvailable(googletest spdlog)

enable_testing()

# Ensure Google Test headers are included
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(
        firstTest
        src/test/ParticleContainerTests.cpp
        src/test/CalculationsTests.cpp
        src/test/ParticleGeneratorTests.cpp
        src/Particle.cpp
        src/Container/ParticleContainer.cpp
        src/Calculations.cpp
        src/ParticleGenerator.cpp
        src/utils/MaxwellBoltzmannDistribution.h
        src/LinkedCell/Cell.cpp
        src/LinkedCell/Cell.h
        src/test/LCParticleContainerTest.cpp
        src/test/CellTest.cpp
        src/Plots/RTCAss3Task2.cpp
        src/Container/LCParticleContainer.cpp
        src/Container/LCParticleContainer.h
        src/Container/BaseParticleContainer.h
        src/test/ThermostatTests.cpp
        src/Thermostat.cpp
        src/Thermostat.h
        src/utils/EpsilonSigma.h
        src/Membrane/Membrane.h
        src/Membrane/MembranePair.h
        src/Membrane/MembranePair.cpp
        src/Membrane/Membrane.cpp
        src/Plots/Ass5Task5Diffusion.cpp
        src/Plots/Ass5Task5RDF.cpp
        src/Plots/Ass5Task2StrategiesSpeedUp.cpp
        src/test/MembraneTest.cpp
)

target_link_libraries(
        firstTest
        GTest::gtest_main
        spdlog::spdlog
)

include(GoogleTest)
include(CTest)
gtest_discover_tests(firstTest)

# Link Google Test libraries with the main executable
target_link_libraries(MolSim PRIVATE GTest::gtest GTest::gtest_main)

set(MEMORYCHECK_COMMAND valgrind)
