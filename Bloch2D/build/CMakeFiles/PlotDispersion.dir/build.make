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
CMAKE_SOURCE_DIR = /home/adhithyar/Documents/Project/EigenCpp/Bloch2D

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/adhithyar/Documents/Project/EigenCpp/Bloch2D/build

# Include any dependencies generated for this target.
include CMakeFiles/PlotDispersion.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/PlotDispersion.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/PlotDispersion.dir/flags.make

CMakeFiles/PlotDispersion.dir/src/plot.cpp.o: CMakeFiles/PlotDispersion.dir/flags.make
CMakeFiles/PlotDispersion.dir/src/plot.cpp.o: ../src/plot.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adhithyar/Documents/Project/EigenCpp/Bloch2D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/PlotDispersion.dir/src/plot.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PlotDispersion.dir/src/plot.cpp.o -c /home/adhithyar/Documents/Project/EigenCpp/Bloch2D/src/plot.cpp

CMakeFiles/PlotDispersion.dir/src/plot.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PlotDispersion.dir/src/plot.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/adhithyar/Documents/Project/EigenCpp/Bloch2D/src/plot.cpp > CMakeFiles/PlotDispersion.dir/src/plot.cpp.i

CMakeFiles/PlotDispersion.dir/src/plot.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PlotDispersion.dir/src/plot.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/adhithyar/Documents/Project/EigenCpp/Bloch2D/src/plot.cpp -o CMakeFiles/PlotDispersion.dir/src/plot.cpp.s

# Object files for target PlotDispersion
PlotDispersion_OBJECTS = \
"CMakeFiles/PlotDispersion.dir/src/plot.cpp.o"

# External object files for target PlotDispersion
PlotDispersion_EXTERNAL_OBJECTS =

PlotDispersion: CMakeFiles/PlotDispersion.dir/src/plot.cpp.o
PlotDispersion: CMakeFiles/PlotDispersion.dir/build.make
PlotDispersion: /usr/lib/x86_64-linux-gnu/libpython3.8.so
PlotDispersion: CMakeFiles/PlotDispersion.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/adhithyar/Documents/Project/EigenCpp/Bloch2D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable PlotDispersion"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PlotDispersion.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/PlotDispersion.dir/build: PlotDispersion

.PHONY : CMakeFiles/PlotDispersion.dir/build

CMakeFiles/PlotDispersion.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/PlotDispersion.dir/cmake_clean.cmake
.PHONY : CMakeFiles/PlotDispersion.dir/clean

CMakeFiles/PlotDispersion.dir/depend:
	cd /home/adhithyar/Documents/Project/EigenCpp/Bloch2D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/adhithyar/Documents/Project/EigenCpp/Bloch2D /home/adhithyar/Documents/Project/EigenCpp/Bloch2D /home/adhithyar/Documents/Project/EigenCpp/Bloch2D/build /home/adhithyar/Documents/Project/EigenCpp/Bloch2D/build /home/adhithyar/Documents/Project/EigenCpp/Bloch2D/build/CMakeFiles/PlotDispersion.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/PlotDispersion.dir/depend

