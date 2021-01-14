#https://evileg.com/en/post/436/

# check requirements ############
cmake_minimum_required(VERSION "3.0")

#settings
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -DNDEBUG")
#set(CMAKE_CXX_FLAGS "-O3")

# project name ##################
project(solver)

add_executable(
	solver 
	src/_main.cpp
)
################################
add_library(
	AutoDiff
	src/AutoDiffScalar.h
	src/AutoDiffJacobian.h
)

set_target_properties(AutoDiff PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(solver PUBLIC AutoDiff)
################################
add_library(
	PortHamiltonian
	src/PortHamiltonian.cpp
	src/PortHamiltonian.h
)

target_link_libraries(solver PUBLIC PortHamiltonian)
################################
add_library(
	PortController
	src/PortController.cpp
	src/PortController.h
)

target_link_libraries(PortHamiltonian PUBLIC PortController)
################################
add_library(
	Cosserat
	src/Cosserat.cpp
	src/Cosserat.h
)

target_link_libraries(PortHamiltonian PUBLIC Cosserat)
################################
add_library(
	LieAlgebra
	src/LieAlgebra.h
)

set_target_properties(LieAlgebra PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(Cosserat PUBLIC LieAlgebra)
################################
add_library(
	ShapeFunctions
	src/Shapes.cpp
	src/Shapes.h
)
target_link_libraries(Cosserat PUBLIC ShapeFunctions)
################################
add_library(
	Config
	src/ConfigFile.cpp
	src/ConfigFile.h
	src/Chameleon.cpp
	src/Chameleon.h
)
target_link_libraries(Cosserat PUBLIC Config)
################################
