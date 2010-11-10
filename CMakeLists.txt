cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


# NeuralGas
PROJECT(NeuralGas)

# Project root directory
GET_FILENAME_COMPONENT(PROJECT_ROOT .. ABSOLUTE)
SET(NEURALGAS_ROOT ${PROJECT_ROOT}/trunk CACHE INTERNAL "Path prefix for NeuralGas package")

# Look for Boost libraries
SET(CMAKE_MODULE_PATH ${NEURALGAS_ROOT})

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
# LINK_LIBRARIES(NeuralGas
  # pthread rt GL glut GLU ${QT_LIBRARIES} qwt-qt4 Ice IceUtil
# )


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


SET(GNG_HEADERS
	${NEURALGAS_ROOT}/GrowingNeuralGas/GNGAlgorithm.h
	${NEURALGAS_ROOT}/GrowingNeuralGas/GNGGraph.h
	${NEURALGAS_ROOT}/GrowingNeuralGas/GNGModul.h
	${NEURALGAS_ROOT}/GrowingNeuralGas/GNGModulGraph.h
)
INSTALL(FILES ${GNG_HEADERS} DESTINATION include/NeuralGas/GrowingNeuralGas/)

SET(ERRORGNG_HEADERS
	${NEURALGAS_ROOT}/GrowingNeuralGas/ErrorBasedGNGAlgorithm/EBGNGAlgorithm.h
)
INSTALL(FILES ${ERRORGNG_HEADERS} DESTINATION include/NeuralGas/GrowingNeuralGas/ErrorBasedGNGAlgorithm)

SET(MGNG_HEADERS
	${NEURALGAS_ROOT}/GrowingNeuralGas/MergeGrowingNeuralGas/MGNGAlgorithm.h
	${NEURALGAS_ROOT}/GrowingNeuralGas/MergeGrowingNeuralGas/MGNGGraph.h
	${NEURALGAS_ROOT}/GrowingNeuralGas/MergeGrowingNeuralGas/CDNAlgorithm.h
)
INSTALL(FILES ${MGNG_HEADERS} DESTINATION include/NeuralGas/GrowingNeuralGas/MergeGrowingNeuralGas)

SET(GRAPH_HEADERS
	${NEURALGAS_ROOT}/Graphs/Base_Graph.h
	${NEURALGAS_ROOT}/Graphs/DGraph.h
	${NEURALGAS_ROOT}/Graphs/TGraph.h
	${NEURALGAS_ROOT}/Graphs/UGraph.h
	${NEURALGAS_ROOT}/Graphs/Vector.h
	${NEURALGAS_ROOT}/Graphs/Vector_temp.h
	${NEURALGAS_ROOT}/Graphs/Vector_temp_easy.h
)
INSTALL(FILES ${GRAPH_HEADERS} DESTINATION include/NeuralGas/Graphs)

SET(DATAGENERATOR_HEADERS
	${NEURALGAS_ROOT}/DataGenerator/DataGenerator.h
	${NEURALGAS_ROOT}/DataGenerator/MackeyGlass.h
)
INSTALL(FILES ${DATAGENERATOR_HEADERS} DESTINATION include/NeuralGas/DataGenerator)

SET(MAINGNG
        ${NEURALGAS_ROOT}/GrowingNeuralGas/mainGNG.cpp
)

SET(MAINNGSUITE
        ${NEURALGAS_ROOT}/mainNGSuite.cpp
)

# NeuralGas library
# ADD_LIBRARY(${PROJECT_NAME} SHARED ${GNG_SOURCES})
# SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE C)
# INSTALL(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib)

ADD_EXECUTABLE(mainGNG ${MAINGNG})
ADD_EXECUTABLE(mainNGSuite ${MAINNGSUITE})

INSTALL(TARGETS mainGNG RUNTIME DESTINATION bin/NeuralGas)
INSTALL(TARGETS mainNGSuite RUNTIME DESTINATION bin/NeuralGas)
