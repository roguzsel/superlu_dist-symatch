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
CMAKE_SOURCE_DIR = /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface

# Include any dependencies generated for this target.
include CMakeFiles/myfort.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/myfort.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/myfort.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/myfort.dir/flags.make

CMakeFiles/myfort.dir/mysub.f.o: CMakeFiles/myfort.dir/flags.make
CMakeFiles/myfort.dir/mysub.f.o: /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/mysub.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/myfort.dir/mysub.f.o"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/mysub.f -o CMakeFiles/myfort.dir/mysub.f.o

CMakeFiles/myfort.dir/mysub.f.i: cmake_force
	@echo "Preprocessing Fortran source to CMakeFiles/myfort.dir/mysub.f.i"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/mysub.f > CMakeFiles/myfort.dir/mysub.f.i

CMakeFiles/myfort.dir/mysub.f.s: cmake_force
	@echo "Compiling Fortran source to assembly CMakeFiles/myfort.dir/mysub.f.s"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/mysub.f -o CMakeFiles/myfort.dir/mysub.f.s

CMakeFiles/myfort.dir/my_sub.f.o: CMakeFiles/myfort.dir/flags.make
CMakeFiles/myfort.dir/my_sub.f.o: /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/my_sub.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/myfort.dir/my_sub.f.o"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/my_sub.f -o CMakeFiles/myfort.dir/my_sub.f.o

CMakeFiles/myfort.dir/my_sub.f.i: cmake_force
	@echo "Preprocessing Fortran source to CMakeFiles/myfort.dir/my_sub.f.i"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/my_sub.f > CMakeFiles/myfort.dir/my_sub.f.i

CMakeFiles/myfort.dir/my_sub.f.s: cmake_force
	@echo "Compiling Fortran source to assembly CMakeFiles/myfort.dir/my_sub.f.s"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/my_sub.f -o CMakeFiles/myfort.dir/my_sub.f.s

CMakeFiles/myfort.dir/mymodule.f90.o: CMakeFiles/myfort.dir/flags.make
CMakeFiles/myfort.dir/mymodule.f90.o: /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/mymodule.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object CMakeFiles/myfort.dir/mymodule.f90.o"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/mymodule.f90 -o CMakeFiles/myfort.dir/mymodule.f90.o

CMakeFiles/myfort.dir/mymodule.f90.i: cmake_force
	@echo "Preprocessing Fortran source to CMakeFiles/myfort.dir/mymodule.f90.i"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/mymodule.f90 > CMakeFiles/myfort.dir/mymodule.f90.i

CMakeFiles/myfort.dir/mymodule.f90.s: cmake_force
	@echo "Compiling Fortran source to assembly CMakeFiles/myfort.dir/mymodule.f90.s"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/mymodule.f90 -o CMakeFiles/myfort.dir/mymodule.f90.s

CMakeFiles/myfort.dir/my_module.f90.o: CMakeFiles/myfort.dir/flags.make
CMakeFiles/myfort.dir/my_module.f90.o: /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/my_module.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object CMakeFiles/myfort.dir/my_module.f90.o"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/my_module.f90 -o CMakeFiles/myfort.dir/my_module.f90.o

CMakeFiles/myfort.dir/my_module.f90.i: cmake_force
	@echo "Preprocessing Fortran source to CMakeFiles/myfort.dir/my_module.f90.i"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/my_module.f90 > CMakeFiles/myfort.dir/my_module.f90.i

CMakeFiles/myfort.dir/my_module.f90.s: cmake_force
	@echo "Compiling Fortran source to assembly CMakeFiles/myfort.dir/my_module.f90.s"
	/opt/cray/pe/craype/2.7.16/bin/ftn $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface/my_module.f90 -o CMakeFiles/myfort.dir/my_module.f90.s

# Object files for target myfort
myfort_OBJECTS = \
"CMakeFiles/myfort.dir/mysub.f.o" \
"CMakeFiles/myfort.dir/my_sub.f.o" \
"CMakeFiles/myfort.dir/mymodule.f90.o" \
"CMakeFiles/myfort.dir/my_module.f90.o"

# External object files for target myfort
myfort_EXTERNAL_OBJECTS =

libmyfort.a: CMakeFiles/myfort.dir/mysub.f.o
libmyfort.a: CMakeFiles/myfort.dir/my_sub.f.o
libmyfort.a: CMakeFiles/myfort.dir/mymodule.f90.o
libmyfort.a: CMakeFiles/myfort.dir/my_module.f90.o
libmyfort.a: CMakeFiles/myfort.dir/build.make
libmyfort.a: CMakeFiles/myfort.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --progress-dir=/global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking Fortran static library libmyfort.a"
	$(CMAKE_COMMAND) -P CMakeFiles/myfort.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/myfort.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/myfort.dir/build: libmyfort.a
.PHONY : CMakeFiles/myfort.dir/build

CMakeFiles/myfort.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/myfort.dir/cmake_clean.cmake
.PHONY : CMakeFiles/myfort.dir/clean

CMakeFiles/myfort.dir/depend:
	cd /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface /global/common/software/nersc/pm-2021q4/sw/cmake-3.22.0/share/cmake-3.22/Modules/FortranCInterface /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface /global/cfs/cdirs/m2956/nanding/myprojects/multi-GPU/superlu_dist/build_perlmutter_1122/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles/myfort.dir/DependInfo.cmake
.PHONY : CMakeFiles/myfort.dir/depend

