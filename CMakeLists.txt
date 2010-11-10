cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


# SMLearning
PROJECT(NeuralGas)


# Project root directory
GET_FILENAME_COMPONENT(PROJECT_ROOT ../.. ABSOLUTE)
SET(NEURALGAS_ROOT ${PROJECT_ROOT}/trunk CACHE INTERNAL "Path prefix for NeuralGas package")

# Look for Boost libraries
SET(CMAKE_MODULE_PATH ${NEURALGAS_ROOT}/cmake)
INCLUDE(${NEURALGAS_ROOT}/cmake/FindBoost.cmake)


# Set the path where other thing should be relative to
GET_FILENAME_COMPONENT(CMAKE_SOURCE_DIR . ABSOLUTE CACHE INTERNAL "")

# Single output directory for building all NeuralGas libraries
SET(LIBRARY_OUTPUT_PATH ${NEURALGAS_ROOT}/lib/linux CACHE PATH "Single output directory for building all NeuralGas libraries")

# Single output directory for building all NeuralGas executables
SET(EXECUTABLE_OUTPUT_PATH ${NEURALGAS_ROOT}/bin/linux CACHE PATH "Single output directory for building all NeuralGas executables")

# Set the install directory.
SET(CMAKE_INSTALL_PREFIX /usr/local CACHE PATH "Path prefix for NeuralGas installation")

# FIND_PACKAGE( Qt4 REQUIRED )

# SET(QWT_INCLUDE /usr/include/qwt-qt4 CACHE PATH "Path prefix for qwt include")
# SET(QWT_LIBRARY /usr/lib/ CACHE PATH "Path prefix for qwt library")

# INCLUDE( ${QT_USE_FILE} )

# SET(ICE_INCLUDE /usr/include CACHE PATH "Path prefix for ZeroC Ice include")
# SET(ICE_LIBRARY /usr/lib CACHE PATH "Path prefix for ZeroC Ice library")
# SET(ICE_BIN /usr/bin CACHE PATH "Path prefix for ZeroC Ice binaries")

ADD_DEFINITIONS(-DLINUX -fpermissive)
LINK_LIBRARIES(NeuralGas
  # pthread rt GL glut GLU ${QT_LIBRARIES} qwt-qt4 Ice IceUtil
)


INCLUDE_DIRECTORIES(
  ${NEURALGAS_ROOT}
#	${QT_INCLUDE_DIR}
#	${QWT_INCLUDE}
#	${ICE_INCLUDE}
)

LINK_DIRECTORIES(
#	${QWT_LIBRARY}
#	${ICE_LIBRARY}
)



ADD_SUBDIRECTORY(${NEURALGAS_ROOT}/DataGenerator)
ADD_SUBDIRECTORY(${NEURALGAS_ROOT}/Graphs)
ADD_SUBDIRECTORY(${NEURALGAS_ROOT}/GrowingNeuralGas)
ADD_SUBDIRECTORY(${NEURALGAS_ROOT}/GrowingNeuralGas/MergeGrowingNeuralGas)
ADD_SUBDIRECTORY(${NEURALGAS_ROOT}/GrowingNeuralGas/Testing)


SET(MAIN_NGSUITE
  ${NEURALGAS_ROOT}/main_NGSuite.cpp
)

# NeuralGas library
ADD_LIBRARY(${PROJECT_NAME} SHARED ${GNG_SOURCES})
INSTALL(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib)

ADD_EXECUTABLE(mainGNG ${MAIN_GNG})
ADD_EXECUTABLE(mainNGSuite ${MAIN_NGSUITE})

INSTALL(TARGETS mainGNG RUNTIME DESTINATION bin/NeuralGas)
INSTALL(TARGETS main_NGSuite RUNTIME DESTINATION bin/NeuralGas)
