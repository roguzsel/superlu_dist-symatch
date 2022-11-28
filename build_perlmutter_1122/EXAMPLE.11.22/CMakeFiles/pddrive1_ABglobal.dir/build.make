# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/bin/cmake

# The command to remove a file.
RM = /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122

# Include any dependencies generated for this target.
include EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/compiler_depend.make

# Include the progress variables for this target.
include EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/progress.make

# Include the compile flags for this target's objects.
include EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/flags.make

EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.o: EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/flags.make
EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.o: ../EXAMPLE/pddrive1_ABglobal.c
EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.o: EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.o"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && /opt/cray/pe/craype/2.7.16/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.o -MF CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.o.d -o CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.o -c /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE/pddrive1_ABglobal.c

EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.i"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && /opt/cray/pe/craype/2.7.16/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE/pddrive1_ABglobal.c > CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.i

EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.s"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && /opt/cray/pe/craype/2.7.16/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE/pddrive1_ABglobal.c -o CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.s

# Object files for target pddrive1_ABglobal
pddrive1_ABglobal_OBJECTS = \
"CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.o"

# External object files for target pddrive1_ABglobal
pddrive1_ABglobal_EXTERNAL_OBJECTS =

EXAMPLE/pddrive1_ABglobal: EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/pddrive1_ABglobal.c.o
EXAMPLE/pddrive1_ABglobal: EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/build.make
EXAMPLE/pddrive1_ABglobal: SRC/libsuperlu_dist.so.6.4.0
EXAMPLE/pddrive1_ABglobal: /opt/cray/pe/libsci/22.06.1.3/GNU/9.1/x86_64/lib/libsci_gnu_81_mp.so
EXAMPLE/pddrive1_ABglobal: /global/cfs/cdirs/m2956/nanding/software/parmetis-4.0.3-perlmutter/build/Linux-x86_64/libparmetis/libparmetis.so
EXAMPLE/pddrive1_ABglobal: /global/cfs/cdirs/m2956/nanding/software/parmetis-4.0.3-perlmutter/build/Linux-x86_64/libmetis/libmetis.so
EXAMPLE/pddrive1_ABglobal: /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64/libcudart.so
EXAMPLE/pddrive1_ABglobal: /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/profilers/Nsight_Compute/../../math_libs/11.7/lib64/libcusolver.so
EXAMPLE/pddrive1_ABglobal: /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/profilers/Nsight_Compute/../../math_libs/11.7/lib64/libcublas.so
EXAMPLE/pddrive1_ABglobal: /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/profilers/Nsight_Compute/../../math_libs/11.7/lib64/libcusparse.so
EXAMPLE/pddrive1_ABglobal: EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable pddrive1_ABglobal"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pddrive1_ABglobal.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/build: EXAMPLE/pddrive1_ABglobal
.PHONY : EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/build

EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/clean:
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && $(CMAKE_COMMAND) -P CMakeFiles/pddrive1_ABglobal.dir/cmake_clean.cmake
.PHONY : EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/clean

EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/depend:
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122 /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : EXAMPLE/CMakeFiles/pddrive1_ABglobal.dir/depend

