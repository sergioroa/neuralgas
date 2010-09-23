cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


# NeuralGas
PROJECT(NeuralGas)


# Project root directory
GET_FILENAME_COMPONENT(PROJECT_ROOT ../.. ABSOLUTE)
SET(NEURALGAS_ROOT ${PROJECT_ROOT}/trunk CACHE INTERNAL "Path prefix for NeuralGas package")
SET(CMAKE_MODULE_PATH ${SMLEARNING_ROOT}/cmake)


# Set the path where other thing should be relative to
GET_FILENAME_COMPONENT(CMAKE_SOURCE_DIR . ABSOLUTE CACHE INTERNAL "")

# Single output directory for building all NeuralGas libraries
SET(LIBRARY_OUTPUT_PATH ${NEURALGAS_ROOT}/output/lib/linux CACHE PATH "Single output directory for building all NeuralGas libraries")

# Single output directory for building all NeuralGas executables
SET(EXECUTABLE_OUTPUT_PATH ${NEURALGAS_ROOT}/output/bin/linux CACHE PATH "Single output directory for building all NeuralGas executables")

# Set the install directory.
SET(CMAKE_INSTALL_PREFIX /usr/local CACHE PATH "Path prefix for NeuralGas installation")

SET(NEURALGAS_INCLUDE ${NEURALGAS_ROOT}/include CACHE PATH "Path prefix for NeuralGas include")

SET(HEADERS
	${NEURALGAS_ROOT}/include/NeuralGas/Base_Graph.h
	${NEURALGAS_ROOT}/include/NeuralGas/DGraph.h
	${NEURALGAS_ROOT}/include/NeuralGas/NeuralGas.h
	${NEURALGAS_ROOT}/include/NeuralGas/TGraph.h
	${NEURALGAS_ROOT}/include/NeuralGas/UGraph.h
	${NEURALGAS_ROOT}/include/NeuralGas/Vector.h
#	${NEURALGAS_ROOT}/include/MGNG/MGNGAlgorithm.h
#	${NEURALGAS_ROOT}/include/MGNG/MGNGGraph.h

)

INSTALL(FILES ${NEURALGAS_HEADERS} DESTINATION include/NeuralGas/)

# NeuralGas library
ADD_LIBRARY(${PROJECT_NAME} SHARED ${HEADERS})
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE C)

INSTALL(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib)
