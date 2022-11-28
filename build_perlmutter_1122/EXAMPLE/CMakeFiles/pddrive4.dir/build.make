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
include EXAMPLE/CMakeFiles/pddrive4.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include EXAMPLE/CMakeFiles/pddrive4.dir/compiler_depend.make

# Include the progress variables for this target.
include EXAMPLE/CMakeFiles/pddrive4.dir/progress.make

# Include the compile flags for this target's objects.
include EXAMPLE/CMakeFiles/pddrive4.dir/flags.make

EXAMPLE/CMakeFiles/pddrive4.dir/pddrive4.c.o: EXAMPLE/CMakeFiles/pddrive4.dir/flags.make
EXAMPLE/CMakeFiles/pddrive4.dir/pddrive4.c.o: ../EXAMPLE/pddrive4.c
EXAMPLE/CMakeFiles/pddrive4.dir/pddrive4.c.o: EXAMPLE/CMakeFiles/pddrive4.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object EXAMPLE/CMakeFiles/pddrive4.dir/pddrive4.c.o"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && /opt/cray/pe/craype/2.7.16/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT EXAMPLE/CMakeFiles/pddrive4.dir/pddrive4.c.o -MF CMakeFiles/pddrive4.dir/pddrive4.c.o.d -o CMakeFiles/pddrive4.dir/pddrive4.c.o -c /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE/pddrive4.c

EXAMPLE/CMakeFiles/pddrive4.dir/pddrive4.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pddrive4.dir/pddrive4.c.i"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && /opt/cray/pe/craype/2.7.16/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE/pddrive4.c > CMakeFiles/pddrive4.dir/pddrive4.c.i

EXAMPLE/CMakeFiles/pddrive4.dir/pddrive4.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pddrive4.dir/pddrive4.c.s"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && /opt/cray/pe/craype/2.7.16/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE/pddrive4.c -o CMakeFiles/pddrive4.dir/pddrive4.c.s

EXAMPLE/CMakeFiles/pddrive4.dir/dcreate_matrix.c.o: EXAMPLE/CMakeFiles/pddrive4.dir/flags.make
EXAMPLE/CMakeFiles/pddrive4.dir/dcreate_matrix.c.o: ../EXAMPLE/dcreate_matrix.c
EXAMPLE/CMakeFiles/pddrive4.dir/dcreate_matrix.c.o: EXAMPLE/CMakeFiles/pddrive4.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object EXAMPLE/CMakeFiles/pddrive4.dir/dcreate_matrix.c.o"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && /opt/cray/pe/craype/2.7.16/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT EXAMPLE/CMakeFiles/pddrive4.dir/dcreate_matrix.c.o -MF CMakeFiles/pddrive4.dir/dcreate_matrix.c.o.d -o CMakeFiles/pddrive4.dir/dcreate_matrix.c.o -c /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE/dcreate_matrix.c

EXAMPLE/CMakeFiles/pddrive4.dir/dcreate_matrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pddrive4.dir/dcreate_matrix.c.i"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && /opt/cray/pe/craype/2.7.16/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE/dcreate_matrix.c > CMakeFiles/pddrive4.dir/dcreate_matrix.c.i

EXAMPLE/CMakeFiles/pddrive4.dir/dcreate_matrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pddrive4.dir/dcreate_matrix.c.s"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && /opt/cray/pe/craype/2.7.16/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE/dcreate_matrix.c -o CMakeFiles/pddrive4.dir/dcreate_matrix.c.s

# Object files for target pddrive4
pddrive4_OBJECTS = \
"CMakeFiles/pddrive4.dir/pddrive4.c.o" \
"CMakeFiles/pddrive4.dir/dcreate_matrix.c.o"

# External object files for target pddrive4
pddrive4_EXTERNAL_OBJECTS =

EXAMPLE/pddrive4: EXAMPLE/CMakeFiles/pddrive4.dir/pddrive4.c.o
EXAMPLE/pddrive4: EXAMPLE/CMakeFiles/pddrive4.dir/dcreate_matrix.c.o
EXAMPLE/pddrive4: EXAMPLE/CMakeFiles/pddrive4.dir/build.make
EXAMPLE/pddrive4: SRC/libsuperlu_dist.so.6.4.0
EXAMPLE/pddrive4: /opt/cray/pe/libsci/22.06.1.3/GNU/9.1/x86_64/lib/libsci_gnu_81_mp.so
EXAMPLE/pddrive4: /global/cfs/cdirs/m2956/nanding/software/parmetis-4.0.3-perlmutter-32bit/build/Linux-x86_64/libparmetis/libparmetis.so
EXAMPLE/pddrive4: /global/cfs/cdirs/m2956/nanding/software/parmetis-4.0.3-perlmutter-32bit/build/Linux-x86_64/libmetis/libmetis.so
EXAMPLE/pddrive4: /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7/lib64/libcudart.so
EXAMPLE/pddrive4: /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/profilers/Nsight_Compute/../../math_libs/11.7/lib64/libcusolver.so
EXAMPLE/pddrive4: /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/profilers/Nsight_Compute/../../math_libs/11.7/lib64/libcublas.so
EXAMPLE/pddrive4: /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/profilers/Nsight_Compute/../../math_libs/11.7/lib64/libcusparse.so
EXAMPLE/pddrive4: EXAMPLE/CMakeFiles/pddrive4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable pddrive4"
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pddrive4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
EXAMPLE/CMakeFiles/pddrive4.dir/build: EXAMPLE/pddrive4
.PHONY : EXAMPLE/CMakeFiles/pddrive4.dir/build

EXAMPLE/CMakeFiles/pddrive4.dir/clean:
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE && $(CMAKE_COMMAND) -P CMakeFiles/pddrive4.dir/cmake_clean.cmake
.PHONY : EXAMPLE/CMakeFiles/pddrive4.dir/clean

EXAMPLE/CMakeFiles/pddrive4.dir/depend:
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/EXAMPLE /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122 /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/EXAMPLE/CMakeFiles/pddrive4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : EXAMPLE/CMakeFiles/pddrive4.dir/depend

