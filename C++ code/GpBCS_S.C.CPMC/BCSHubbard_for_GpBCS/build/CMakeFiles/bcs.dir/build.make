# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build

# Include any dependencies generated for this target.
include CMakeFiles/bcs.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/bcs.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/bcs.dir/flags.make

CMakeFiles/bcs.dir/source/bcs.cpp.o: CMakeFiles/bcs.dir/flags.make
CMakeFiles/bcs.dir/source/bcs.cpp.o: ../source/bcs.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/bcs.dir/source/bcs.cpp.o"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bcs.dir/source/bcs.cpp.o -c /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcs.cpp

CMakeFiles/bcs.dir/source/bcs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bcs.dir/source/bcs.cpp.i"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcs.cpp > CMakeFiles/bcs.dir/source/bcs.cpp.i

CMakeFiles/bcs.dir/source/bcs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bcs.dir/source/bcs.cpp.s"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcs.cpp -o CMakeFiles/bcs.dir/source/bcs.cpp.s

CMakeFiles/bcs.dir/source/bcsDetail.cpp.o: CMakeFiles/bcs.dir/flags.make
CMakeFiles/bcs.dir/source/bcsDetail.cpp.o: ../source/bcsDetail.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/bcs.dir/source/bcsDetail.cpp.o"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bcs.dir/source/bcsDetail.cpp.o -c /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsDetail.cpp

CMakeFiles/bcs.dir/source/bcsDetail.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bcs.dir/source/bcsDetail.cpp.i"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsDetail.cpp > CMakeFiles/bcs.dir/source/bcsDetail.cpp.i

CMakeFiles/bcs.dir/source/bcsDetail.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bcs.dir/source/bcsDetail.cpp.s"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsDetail.cpp -o CMakeFiles/bcs.dir/source/bcsDetail.cpp.s

CMakeFiles/bcs.dir/source/bcsMethod.cpp.o: CMakeFiles/bcs.dir/flags.make
CMakeFiles/bcs.dir/source/bcsMethod.cpp.o: ../source/bcsMethod.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/bcs.dir/source/bcsMethod.cpp.o"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bcs.dir/source/bcsMethod.cpp.o -c /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsMethod.cpp

CMakeFiles/bcs.dir/source/bcsMethod.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bcs.dir/source/bcsMethod.cpp.i"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsMethod.cpp > CMakeFiles/bcs.dir/source/bcsMethod.cpp.i

CMakeFiles/bcs.dir/source/bcsMethod.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bcs.dir/source/bcsMethod.cpp.s"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsMethod.cpp -o CMakeFiles/bcs.dir/source/bcsMethod.cpp.s

CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.o: CMakeFiles/bcs.dir/flags.make
CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.o: ../source/bcsTuneMu.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.o"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.o -c /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsTuneMu.cpp

CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.i"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsTuneMu.cpp > CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.i

CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.s"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsTuneMu.cpp -o CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.s

CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.o: CMakeFiles/bcs.dir/flags.make
CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.o: ../source/bcsTunePairing.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.o"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.o -c /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsTunePairing.cpp

CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.i"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsTunePairing.cpp > CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.i

CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.s"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/bcsTunePairing.cpp -o CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.s

CMakeFiles/bcs.dir/source/main.cpp.o: CMakeFiles/bcs.dir/flags.make
CMakeFiles/bcs.dir/source/main.cpp.o: ../source/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/bcs.dir/source/main.cpp.o"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bcs.dir/source/main.cpp.o -c /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/main.cpp

CMakeFiles/bcs.dir/source/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bcs.dir/source/main.cpp.i"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/main.cpp > CMakeFiles/bcs.dir/source/main.cpp.i

CMakeFiles/bcs.dir/source/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bcs.dir/source/main.cpp.s"
	/home/xiaozhiyu/Desktop/personal_lib/mpi/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/source/main.cpp -o CMakeFiles/bcs.dir/source/main.cpp.s

# Object files for target bcs
bcs_OBJECTS = \
"CMakeFiles/bcs.dir/source/bcs.cpp.o" \
"CMakeFiles/bcs.dir/source/bcsDetail.cpp.o" \
"CMakeFiles/bcs.dir/source/bcsMethod.cpp.o" \
"CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.o" \
"CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.o" \
"CMakeFiles/bcs.dir/source/main.cpp.o"

# External object files for target bcs
bcs_EXTERNAL_OBJECTS =

bcs: CMakeFiles/bcs.dir/source/bcs.cpp.o
bcs: CMakeFiles/bcs.dir/source/bcsDetail.cpp.o
bcs: CMakeFiles/bcs.dir/source/bcsMethod.cpp.o
bcs: CMakeFiles/bcs.dir/source/bcsTuneMu.cpp.o
bcs: CMakeFiles/bcs.dir/source/bcsTunePairing.cpp.o
bcs: CMakeFiles/bcs.dir/source/main.cpp.o
bcs: CMakeFiles/bcs.dir/build.make
bcs: /home/xiaozhiyu/lib/afqmclab/lib/libafqmcHao.a
bcs: /home/xiaozhiyu/lib/afqmclab/lib/liblanczosHao.a
bcs: /home/xiaozhiyu/lib/afqmclab/lib/libcommonHao.a
bcs: /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.so
bcs: /opt/intel/mkl/lib/intel64/libmkl_core.so
bcs: /opt/intel/mkl/lib/intel64/libmkl_sequential.so
bcs: /usr/local/lib/libfftw3.a
bcs: /home/xiaozhiyu/Desktop/personal_lib/sprng2.0/lib/libsprng.a
bcs: /usr/local/lib/libgmp.a
bcs: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.so
bcs: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so
bcs: /usr/lib/x86_64-linux-gnu/libpthread.so
bcs: /usr/lib/x86_64-linux-gnu/libsz.so
bcs: /usr/lib/x86_64-linux-gnu/libz.so
bcs: /usr/lib/x86_64-linux-gnu/libdl.so
bcs: /usr/lib/x86_64-linux-gnu/libm.so
bcs: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.so
bcs: /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.so
bcs: /usr/local/lib/libgtest.so
bcs: CMakeFiles/bcs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable bcs"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bcs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/bcs.dir/build: bcs

.PHONY : CMakeFiles/bcs.dir/build

CMakeFiles/bcs.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/bcs.dir/cmake_clean.cmake
.PHONY : CMakeFiles/bcs.dir/clean

CMakeFiles/bcs.dir/depend:
	cd /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build /home/xiaozhiyu/Desktop/code/code/C++/AFQMC_XIao_from_Hao_code/BCSHubbard/build/CMakeFiles/bcs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/bcs.dir/depend
