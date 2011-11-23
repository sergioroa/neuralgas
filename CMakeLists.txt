cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


# NeuralGas
PROJECT(NeuralGas)

# Project root directory
GET_FILENAME_COMPONENT(PROJECT_ROOT .. ABSOLUTE)
SET(NEURALGAS_ROOT ${PROJECT_ROOT}/trunk CACHE INTERNAL "Path prefix for NeuralGas package")

## Look for Boost libraries
# SET(CMAKE_MODULE_PATH ${NEURALGAS_ROOT})
# INCLUDE(${NEURALGAS_ROOT}/FindBoost.cmake)

# Set the path where other thing should be relative to
GET_FILENAME_COMPONENT(CMAKE_SOURCE_DIR . ABSOLUTE CACHE INTERNAL "")

SET(Boost_MAJOR_VERSION 1)
SET(Boost_MINOR_VERSION 39)
FIND_PACKAGE( Boost 1.39 COMPONENTS program_options )

# Single output directory for building all NeuralGas libraries
SET(LIBRARY_OUTPUT_PATH ${NEURALGAS_ROOT}/lib/linux CACHE PATH "Single output directory for building all NeuralGas libraries")

# Single output directory for building all NeuralGas executables
SET(EXECUTABLE_OUTPUT_PATH ${NEURALGAS_ROOT}/bin/linux CACHE PATH "Single output directory for building all NeuralGas executables")

# Set the install directory.
SET(CMAKE_INSTALL_PREFIX /usr/local CACHE PATH "Path prefix for NeuralGas installation")

FIND_PACKAGE( Qt4 REQUIRED )

# SET(QWT_INCLUDE /usr/include/qwt-qt4 CACHE PATH "Path prefix for qwt include")
# SET(QWT_LIBRARY /usr/lib/ CACHE PATH "Path prefix for qwt library")

INCLUDE( ${QT_USE_FILE} )

# SET(ICE_INCLUDE /usr/include CACHE PATH "Path prefix for ZeroC Ice include")
# SET(ICE_LIBRARY /usr/lib CACHE PATH "Path prefix for ZeroC Ice library")
# SET(ICE_BIN /usr/bin CACHE PATH "Path prefix for ZeroC Ice binaries")

# ADD_DEFINITIONS(-DLINUX -fpermissive)
ADD_DEFINITIONS(-Wall -DLINUX)

LINK_LIBRARIES(NeuralGas ${QT_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY}
  # pthread rt GL glut GLU qwt-qt4 Ice IceUtil
)


INCLUDE_DIRECTORIES(
        ${NEURALGAS_ROOT}
	${QT_INCLUDE_DIR}
	${Boost_INCLUDE_DIRS}
#	${QWT_INCLUDE}
#	${ICE_INCLUDE}
)

LINK_DIRECTORIES(
	${Boost_LIBRARY_DIRS}
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
INSTALL(FILES ${ERRORGNG_HEADERS} DESTINATION include/NeuralGas/GrowingNeuralGas/ErrorBasedGNGAlgorithm/)

SET(LLRGNG_HEADERS
	${NEURALGAS_ROOT}/GrowingNeuralGas/LifelongRobustGNGAlgorithm/LLRGNGAlgorithm.h
	${NEURALGAS_ROOT}/GrowingNeuralGas/LifelongRobustGNGAlgorithm/LLRGNGGraph.h
)
INSTALL(FILES ${LLRGNG_HEADERS} DESTINATION include/NeuralGas/GrowingNeuralGas/LifelongRobustGNGAlgorithm/)

SET(LLRGNG_MOC_HEADERS
 	${NEURALGAS_ROOT}/GrowingNeuralGas/LifelongRobustGNGAlgorithm/LLRGNGAlgorithm.h
)
QT4_WRAP_CPP(LLRGNG_GEN_MOC ${LLRGNG_MOC_HEADERS} )

SET(LLRGNG_SOURCES
        ${LLRGNG_GEN_MOC}
)

SET(MGNG_HEADERS
	${NEURALGAS_ROOT}/GrowingNeuralGas/MergeGrowingNeuralGas/MGNGAlgorithm.h
	${NEURALGAS_ROOT}/GrowingNeuralGas/MergeGrowingNeuralGas/MGNGGraph.h
	${NEURALGAS_ROOT}/GrowingNeuralGas/MergeGrowingNeuralGas/CDNAlgorithm.h
)
INSTALL(FILES ${MGNG_HEADERS} DESTINATION include/NeuralGas/GrowingNeuralGas/MergeGrowingNeuralGas/)

SET(GNGTESTING_HEADERS
	${NEURALGAS_ROOT}/GrowingNeuralGas/Testing/ErrorTesting.h)
INSTALL(FILES ${MGNG_HEADERS} DESTINATION include/NeuralGas/GrowingNeuralGas/Testing/)

SET(GRAPH_HEADERS
	${NEURALGAS_ROOT}/Graphs/Base_Graph.h
	${NEURALGAS_ROOT}/Graphs/DGraph.h
	${NEURALGAS_ROOT}/Graphs/TGraph.h
	${NEURALGAS_ROOT}/Graphs/UGraph.h
	${NEURALGAS_ROOT}/Graphs/Vector.h
)
INSTALL(FILES ${GRAPH_HEADERS} DESTINATION include/NeuralGas/Graphs/)

SET(DATAGENERATOR_HEADERS
	${NEURALGAS_ROOT}/DataGenerator/DataGenerator.h
	${NEURALGAS_ROOT}/DataGenerator/MackeyGlass.h
	${NEURALGAS_ROOT}/DataGenerator/BinaryAutomata.h
	${NEURALGAS_ROOT}/DataGenerator/NoisyAutomata.h
	${NEURALGAS_ROOT}/DataGenerator/GaussianNoise.h
	${NEURALGAS_ROOT}/DataGenerator/ProbMackeyGlass.h
)
INSTALL(FILES ${DATAGENERATOR_HEADERS} DESTINATION include/NeuralGas/DataGenerator/)

SET(VORONOIDIAGRAMM_HEADERS
	${NEURALGAS_ROOT}/VoronoiDiagramm/VoronoiDiagramGenerator.h
	${NEURALGAS_ROOT}/VoronoiDiagramm/Voronoi.h
	${NEURALGAS_ROOT}/VoronoiDiagramm/VoronoiWidget.h)
INSTALL(FILES ${VORONOIDIAGRAMM_HEADERS} DESTINATION include/NeuralGas/VoronoiDiagramm/)

SET(TOOLS_HEADERS
	${NEURALGAS_ROOT}/tools/metrics.h
	${NEURALGAS_ROOT}/tools/math_helpers.h
	${NEURALGAS_ROOT}/tools/helpers.h)
INSTALL(FILES ${TOOLS_HEADERS} DESTINATION include/NeuralGas/tools/)

SET(VORONOIDIAGRAMM_MOC_HEADERS
	${NEURALGAS_ROOT}/VoronoiDiagramm/VoronoiWidget.h
)
QT4_WRAP_CPP(VORONOIDIAGRAMM_GEN_MOC ${VORONOIDIAGRAMM_MOC_HEADERS} )

SET(VORONOIDIAGRAMM_SOURCES
	${NEURALGAS_ROOT}/VoronoiDiagramm/VoronoiDiagramGenerator.cpp
	${NEURALGAS_ROOT}/VoronoiDiagramm/Voronoi.cpp
	${NEURALGAS_ROOT}/VoronoiDiagramm/VoronoiWidget.cpp
	${VORONOIDIAGRAMM_GEN_MOC}
)

SET(HEADERS
	${NEURALGAS_ROOT}/NeuralGas.h
	${NEURALGAS_ROOT}/NeuralGasSuite.h)
INSTALL(FILES ${HEADERS} DESTINATION include/NeuralGas/)

SET(MACKEYGLASS
        ${NEURALGAS_ROOT}/DataGenerator/mainMackeyGlass.cpp)

SET(BINARYAUTOMATA
        ${NEURALGAS_ROOT}/DataGenerator/mainBinaryAutomata.cpp)

SET(PROBMACKEYGLASS
        ${NEURALGAS_ROOT}/DataGenerator/mainProbMackeyGlass.cpp)

SET(NOISYAUTOMATA
        ${NEURALGAS_ROOT}/DataGenerator/mainNoisyAutomata.cpp)

SET(GAUSSIANNOISE
        ${NEURALGAS_ROOT}/DataGenerator/mainGaussianNoise.cpp)

SET(VORONOIDIAGRAMM
        ${NEURALGAS_ROOT}/VoronoiDiagramm/main.cpp)

SET(VORONOIWIDGET
        ${NEURALGAS_ROOT}/VoronoiDiagramm/mainVoronoiWidget.cpp)

SET(LLRGNG
        ${NEURALGAS_ROOT}/GrowingNeuralGas/LifelongRobustGNGAlgorithm/mainLLRGNG.cpp)

SET(EBGNG
        ${NEURALGAS_ROOT}/GrowingNeuralGas/ErrorBasedGNGAlgorithm/mainEBGNG.cpp)

SET(MGNG
        ${NEURALGAS_ROOT}/GrowingNeuralGas/MergeGrowingNeuralGas/MGNGmain.cpp)

SET(CDN
        ${NEURALGAS_ROOT}/GrowingNeuralGas/MergeGrowingNeuralGas/CDNmain.cpp)

SET(GNG
        ${NEURALGAS_ROOT}/GrowingNeuralGas/mainGNG.cpp)

SET(GNGTESTING
        ${NEURALGAS_ROOT}/GrowingNeuralGas/Testing/mainTesting.cpp)

SET(GRAPHS
        ${NEURALGAS_ROOT}/Graphs/main.cpp)

SET(NGSUITE
        ${NEURALGAS_ROOT}/mainNGSuite.cpp)

# NeuralGas library
# ADD_LIBRARY(${PROJECT_NAME} SHARED ${GNG_SOURCES})
ADD_LIBRARY(${PROJECT_NAME} SHARED ${VORONOIDIAGRAMM_SOURCES}
${LLRGNG_SOURCES}
)
# SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE C)
INSTALL(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib)

ADD_EXECUTABLE(NGSuite ${NGSUITE})
ADD_EXECUTABLE(MackeyGlass ${MACKEYGLASS})
ADD_EXECUTABLE(BinaryAutomata ${BINARYAUTOMATA})
ADD_EXECUTABLE(ProbMackeyGlass ${PROBMACKEYGLASS})
ADD_EXECUTABLE(NoisyAutomata ${NOISYAUTOMATA})
ADD_EXECUTABLE(GaussianNoise ${GAUSSIANNOISE})
ADD_EXECUTABLE(VoronoiDiagramm ${VORONOIDIAGRAMM})
ADD_EXECUTABLE(VoronoiWidget ${VORONOIWIDGET})
ADD_EXECUTABLE(GNG ${GNG})
ADD_EXECUTABLE(LLRGNG ${LLRGNG})
ADD_EXECUTABLE(EBGNG ${EBGNG})
ADD_EXECUTABLE(MGNG ${MGNG})
ADD_EXECUTABLE(CDN ${CDN})
ADD_EXECUTABLE(GNGTesting ${GNGTESTING})
ADD_EXECUTABLE(graphs ${GRAPHS})


# INSTALL(TARGETS GNG RUNTIME DESTINATION bin/NeuralGas)
INSTALL(TARGETS NGSuite RUNTIME DESTINATION bin/NeuralGas)
INSTALL(TARGETS LLRGNG RUNTIME DESTINATION bin/NeuralGas)
INSTALL(TARGETS GaussianNoise RUNTIME DESTINATION bin/NeuralGas)
INSTALL(TARGETS EBGNG RUNTIME DESTINATION bin/NeuralGas)
INSTALL(TARGETS GNGTesting RUNTIME DESTINATION bin/NeuralGas)
INSTALL(TARGETS VoronoiDiagramm RUNTIME DESTINATION bin/NeuralGas)
INSTALL(TARGETS VoronoiWidget RUNTIME DESTINATION bin/NeuralGas)
