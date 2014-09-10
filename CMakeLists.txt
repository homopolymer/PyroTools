# ================================================ +
# PyroTools -- A toolkit of computational programs +
#              for sequencing data analysis        +
#                                                  +
# Top-level                                        +
# ================================================ +

# set project name
project(PyroTools)

# cmake requirement
cmake_minimum_required(VERSION 2.8)

# set binary directory
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")

# ================================================ +
# External Projects                                +
# ================================================ +
include(ExternalProject)

# External Project -- BamTools
#find_package(ZLIB REQUIRED)
set(BamTools_DIR ${CMAKE_SOURCE_DIR}/bamtools)
set(BamTools_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${BamTools_DIR})
ExternalProject_Add(BamTools 
    PREFIX ${BamTools_DIR}
    DOWNLOAD_COMMAND ""
    UPDATE_COMMAND ""
    DOWNLOAD_DIR ${BamTools_DIR}
    SOURCE_DIR ${BamTools_DIR}
    BINARY_DIR ${BamTools_DIR}/build
    INSTALL_DIR ${BamTools_DIR}
    CMAKE_ARGS ${BamTools_CMAKE_ARGS}
)

# External Project -- NLOPT
set(Nlopt_DIR ${CMAKE_SOURCE_DIR}/nlopt)
ExternalProject_Add(Nlopt
    PREFIX ${Nlopt_DIR}
    DOWNLOAD_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ${Nlopt_DIR}/configure --prefix=${Nlopt_DIR}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)

# ================================================ +
# Compile setting for PyroTools                    +
# ================================================ +

# define compiler flags for all code
set(CMAKE_BUILD_TYPE Debug)
add_definitions(-O3 -Wall -D_FILE_OFFSET_BITS=64 -std=c++11)

# add our include root path
include_directories(src)

# list subdirectory to build in
add_subdirectory(src)

#aux_source_directory(. SRC_LIST)
#add_executable(${PROJECT_NAME} ${SRC_LIST})


