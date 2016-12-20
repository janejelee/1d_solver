# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs

# Include any dependencies generated for this target.
include CMakeFiles/dofs.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/dofs.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dofs.dir/flags.make

CMakeFiles/dofs.dir/dofs.cc.o: CMakeFiles/dofs.dir/flags.make
CMakeFiles/dofs.dir/dofs.cc.o: dofs.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/dofs.dir/dofs.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dofs.dir/dofs.cc.o -c /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs/dofs.cc

CMakeFiles/dofs.dir/dofs.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dofs.dir/dofs.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs/dofs.cc > CMakeFiles/dofs.dir/dofs.cc.i

CMakeFiles/dofs.dir/dofs.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dofs.dir/dofs.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs/dofs.cc -o CMakeFiles/dofs.dir/dofs.cc.s

CMakeFiles/dofs.dir/dofs.cc.o.requires:

.PHONY : CMakeFiles/dofs.dir/dofs.cc.o.requires

CMakeFiles/dofs.dir/dofs.cc.o.provides: CMakeFiles/dofs.dir/dofs.cc.o.requires
	$(MAKE) -f CMakeFiles/dofs.dir/build.make CMakeFiles/dofs.dir/dofs.cc.o.provides.build
.PHONY : CMakeFiles/dofs.dir/dofs.cc.o.provides

CMakeFiles/dofs.dir/dofs.cc.o.provides.build: CMakeFiles/dofs.dir/dofs.cc.o


# Object files for target dofs
dofs_OBJECTS = \
"CMakeFiles/dofs.dir/dofs.cc.o"

# External object files for target dofs
dofs_EXTERNAL_OBJECTS =

dofs: CMakeFiles/dofs.dir/dofs.cc.o
dofs: CMakeFiles/dofs.dir/build.make
dofs: /scratch/leej/deal.II/lib/libdeal_II.g.so.8.4.1
dofs: /usr/lib/x86_64-linux-gnu/libbz2.so
dofs: /usr/lib/openmpi/lib/libmpi_usempif08.so
dofs: /usr/lib/openmpi/lib/libmpi_usempi_ignore_tkr.so
dofs: /usr/lib/openmpi/lib/libmpi_mpifh.so
dofs: /usr/lib/x86_64-linux-gnu/libz.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libmuelu-adapters.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libmuelu-interface.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libmuelu.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libteko.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libstratimikos.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libstratimikosbelos.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libstratimikosaztecoo.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libstratimikosamesos.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libstratimikosml.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libstratimikosifpack.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libifpack2-adapters.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libifpack2.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libanasazitpetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libModeLaplace.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libanasaziepetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libanasazi.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libamesos2.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libbelostpetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libbelosepetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libbelos.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libml.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libifpack.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libzoltan2.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libpamgen_extras.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libpamgen.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libamesos.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libgaleri-xpetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libgaleri-epetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libaztecoo.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libisorropia.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libxpetra-sup.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libxpetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libthyratpetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libthyraepetraext.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libthyraepetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libthyracore.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libepetraext.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libtpetraext.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libtpetrainout.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libtpetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libkokkostsqr.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libtpetrakernels.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libtpetraclassiclinalg.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libtpetraclassicnodeapi.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libtpetraclassic.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libtriutils.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libzoltan.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libepetra.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libsacado.so
dofs: /scratch/trilinos-12.10.1-Source/lib/librtop.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libteuchoskokkoscomm.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libteuchoskokkoscompat.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libteuchosremainder.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libteuchosnumerics.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libteuchoscomm.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libteuchosparameterlist.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libteuchoscore.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libkokkosalgorithms.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libkokkoscontainers.so
dofs: /scratch/trilinos-12.10.1-Source/lib/libkokkoscore.so
dofs: /usr/lib/libblas.so
dofs: /usr/lib/openmpi/lib/libmpi_cxx.so
dofs: /usr/lib/x86_64-linux-gnu/libumfpack.so
dofs: /usr/lib/x86_64-linux-gnu/libcholmod.so
dofs: /usr/lib/x86_64-linux-gnu/libccolamd.so
dofs: /usr/lib/x86_64-linux-gnu/libcolamd.so
dofs: /usr/lib/x86_64-linux-gnu/libcamd.so
dofs: /usr/lib/x86_64-linux-gnu/libamd.so
dofs: /usr/lib/x86_64-linux-gnu/libmetis.so
dofs: /usr/lib/libparpack.so
dofs: /usr/lib/libarpack.so
dofs: /usr/lib/liblapack.so
dofs: /usr/lib/libf77blas.so
dofs: /usr/lib/libatlas.so
dofs: /usr/lib/openmpi/lib/libmpi.so
dofs: /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so
dofs: /usr/lib/x86_64-linux-gnu/libnetcdf.so
dofs: /usr/lib/x86_64-linux-gnu/libslepc.so
dofs: /scratch/leej/petsc-3.6.4/x86_64/lib/libpetsc.so
dofs: CMakeFiles/dofs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable dofs"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dofs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dofs.dir/build: dofs

.PHONY : CMakeFiles/dofs.dir/build

CMakeFiles/dofs.dir/requires: CMakeFiles/dofs.dir/dofs.cc.o.requires

.PHONY : CMakeFiles/dofs.dir/requires

CMakeFiles/dofs.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dofs.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dofs.dir/clean

CMakeFiles/dofs.dir/depend:
	cd /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs /scratch/leej/deal.II/my_examples/my_step-33/tests/test_dofs/CMakeFiles/dofs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dofs.dir/depend

