# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/build_continious

# Include any dependencies generated for this target.
include CMakeFiles/solver.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/solver.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/solver.dir/flags.make

CMakeFiles/solver.dir/src/_main.cpp.o: CMakeFiles/solver.dir/flags.make
CMakeFiles/solver.dir/src/_main.cpp.o: ../src/_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/build_continious/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/solver.dir/src/_main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/solver.dir/src/_main.cpp.o -c /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/src/_main.cpp

CMakeFiles/solver.dir/src/_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/solver.dir/src/_main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/src/_main.cpp > CMakeFiles/solver.dir/src/_main.cpp.i

CMakeFiles/solver.dir/src/_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/solver.dir/src/_main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/src/_main.cpp -o CMakeFiles/solver.dir/src/_main.cpp.s

# Object files for target solver
solver_OBJECTS = \
"CMakeFiles/solver.dir/src/_main.cpp.o"

# External object files for target solver
solver_EXTERNAL_OBJECTS =

solver: CMakeFiles/solver.dir/src/_main.cpp.o
solver: CMakeFiles/solver.dir/build.make
solver: libAutoDiff.a
solver: libPortHamiltonian.a
solver: libPortController.a
solver: libCosserat.a
solver: libLieAlgebra.a
solver: libShapeFunctions.a
solver: libConfig.a
solver: CMakeFiles/solver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/build_continious/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable solver"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/solver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/solver.dir/build: solver

.PHONY : CMakeFiles/solver.dir/build

CMakeFiles/solver.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/solver.dir/cmake_clean.cmake
.PHONY : CMakeFiles/solver.dir/clean

CMakeFiles/solver.dir/depend:
	cd /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/build_continious && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/build_continious /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/build_continious /home/brandon/Documents/MATLAB/SorotokiCode/src/model/tools/solver/build_continious/CMakeFiles/solver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/solver.dir/depend

