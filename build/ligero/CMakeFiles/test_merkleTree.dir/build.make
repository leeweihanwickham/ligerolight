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
CMAKE_SOURCE_DIR = "/home/leeweihan/Desktop/ligero++ -1209"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/leeweihan/Desktop/ligero++ -1209/build"

# Include any dependencies generated for this target.
include ligero/CMakeFiles/test_merkleTree.dir/depend.make

# Include the progress variables for this target.
include ligero/CMakeFiles/test_merkleTree.dir/progress.make

# Include the compile flags for this target's objects.
include ligero/CMakeFiles/test_merkleTree.dir/flags.make

ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o: ligero/CMakeFiles/test_merkleTree.dir/flags.make
ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o: ../ligero/tests/test_merkleTree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/leeweihan/Desktop/ligero++ -1209/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o -c "/home/leeweihan/Desktop/ligero++ -1209/ligero/tests/test_merkleTree.cpp"

ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.i"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/leeweihan/Desktop/ligero++ -1209/ligero/tests/test_merkleTree.cpp" > CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.i

ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.s"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/leeweihan/Desktop/ligero++ -1209/ligero/tests/test_merkleTree.cpp" -o CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.s

ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o.requires:

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o.requires

ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o.provides: ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o.requires
	$(MAKE) -f ligero/CMakeFiles/test_merkleTree.dir/build.make ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o.provides.build
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o.provides

ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o.provides.build: ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o


ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o: ligero/CMakeFiles/test_merkleTree.dir/flags.make
ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o: ../ligero/bcs/BLAKE3/blake3.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/leeweihan/Desktop/ligero++ -1209/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building C object ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o   -c "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3.c"

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.i"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3.c" > CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.i

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.s"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3.c" -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.s

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o.requires:

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o.requires

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o.provides: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o.requires
	$(MAKE) -f ligero/CMakeFiles/test_merkleTree.dir/build.make ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o.provides.build
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o.provides

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o.provides.build: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o


ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o: ligero/CMakeFiles/test_merkleTree.dir/flags.make
ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o: ../ligero/bcs/BLAKE3/blake3_dispatch.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/leeweihan/Desktop/ligero++ -1209/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building C object ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o   -c "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_dispatch.c"

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.i"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_dispatch.c" > CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.i

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.s"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_dispatch.c" -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.s

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o.requires:

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o.requires

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o.provides: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o.requires
	$(MAKE) -f ligero/CMakeFiles/test_merkleTree.dir/build.make ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o.provides.build
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o.provides

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o.provides.build: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o


ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o: ligero/CMakeFiles/test_merkleTree.dir/flags.make
ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o: ../ligero/bcs/BLAKE3/blake3_portable.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/leeweihan/Desktop/ligero++ -1209/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building C object ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o   -c "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_portable.c"

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.i"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_portable.c" > CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.i

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.s"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_portable.c" -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.s

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o.requires:

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o.requires

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o.provides: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o.requires
	$(MAKE) -f ligero/CMakeFiles/test_merkleTree.dir/build.make ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o.provides.build
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o.provides

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o.provides.build: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o


ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o: ligero/CMakeFiles/test_merkleTree.dir/flags.make
ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o: ../ligero/bcs/BLAKE3/blake3_sse2_x86-64_unix.S
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/leeweihan/Desktop/ligero++ -1209/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Building ASM object ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(ASM_DEFINES) $(ASM_INCLUDES) $(ASM_FLAGS) -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o -c "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_sse2_x86-64_unix.S"

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o.requires:

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o.requires

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o.provides: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o.requires
	$(MAKE) -f ligero/CMakeFiles/test_merkleTree.dir/build.make ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o.provides.build
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o.provides

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o.provides.build: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o


ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o: ligero/CMakeFiles/test_merkleTree.dir/flags.make
ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o: ../ligero/bcs/BLAKE3/blake3_sse41_x86-64_unix.S
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/leeweihan/Desktop/ligero++ -1209/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Building ASM object ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(ASM_DEFINES) $(ASM_INCLUDES) $(ASM_FLAGS) -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o -c "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_sse41_x86-64_unix.S"

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o.requires:

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o.requires

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o.provides: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o.requires
	$(MAKE) -f ligero/CMakeFiles/test_merkleTree.dir/build.make ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o.provides.build
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o.provides

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o.provides.build: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o


ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o: ligero/CMakeFiles/test_merkleTree.dir/flags.make
ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o: ../ligero/bcs/BLAKE3/blake3_avx2_x86-64_unix.S
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/leeweihan/Desktop/ligero++ -1209/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_7) "Building ASM object ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(ASM_DEFINES) $(ASM_INCLUDES) $(ASM_FLAGS) -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o -c "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_avx2_x86-64_unix.S"

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o.requires:

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o.requires

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o.provides: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o.requires
	$(MAKE) -f ligero/CMakeFiles/test_merkleTree.dir/build.make ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o.provides.build
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o.provides

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o.provides.build: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o


ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o: ligero/CMakeFiles/test_merkleTree.dir/flags.make
ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o: ../ligero/bcs/BLAKE3/blake3_avx512_x86-64_unix.S
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/leeweihan/Desktop/ligero++ -1209/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_8) "Building ASM object ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && /usr/bin/cc  $(ASM_DEFINES) $(ASM_INCLUDES) $(ASM_FLAGS) -o CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o -c "/home/leeweihan/Desktop/ligero++ -1209/ligero/bcs/BLAKE3/blake3_avx512_x86-64_unix.S"

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o.requires:

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o.requires

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o.provides: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o.requires
	$(MAKE) -f ligero/CMakeFiles/test_merkleTree.dir/build.make ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o.provides.build
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o.provides

ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o.provides.build: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o


# Object files for target test_merkleTree
test_merkleTree_OBJECTS = \
"CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o" \
"CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o" \
"CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o" \
"CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o" \
"CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o" \
"CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o" \
"CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o" \
"CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o"

# External object files for target test_merkleTree
test_merkleTree_EXTERNAL_OBJECTS =

ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o
ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o
ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o
ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o
ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o
ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o
ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o
ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o
ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/build.make
ligero/test_merkleTree: ligero/libligero.a
ligero/test_merkleTree: depends/libff/libff/libff.a
ligero/test_merkleTree: /usr/lib/x86_64-linux-gnu/libgmp.so
ligero/test_merkleTree: depends/libzm.a
ligero/test_merkleTree: /usr/lib/x86_64-linux-gnu/libsodium.so
ligero/test_merkleTree: ligero/CMakeFiles/test_merkleTree.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/leeweihan/Desktop/ligero++ -1209/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable test_merkleTree"
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_merkleTree.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ligero/CMakeFiles/test_merkleTree.dir/build: ligero/test_merkleTree

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/build

ligero/CMakeFiles/test_merkleTree.dir/requires: ligero/CMakeFiles/test_merkleTree.dir/tests/test_merkleTree.cpp.o.requires
ligero/CMakeFiles/test_merkleTree.dir/requires: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3.c.o.requires
ligero/CMakeFiles/test_merkleTree.dir/requires: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_dispatch.c.o.requires
ligero/CMakeFiles/test_merkleTree.dir/requires: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_portable.c.o.requires
ligero/CMakeFiles/test_merkleTree.dir/requires: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse2_x86-64_unix.S.o.requires
ligero/CMakeFiles/test_merkleTree.dir/requires: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_sse41_x86-64_unix.S.o.requires
ligero/CMakeFiles/test_merkleTree.dir/requires: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx2_x86-64_unix.S.o.requires
ligero/CMakeFiles/test_merkleTree.dir/requires: ligero/CMakeFiles/test_merkleTree.dir/bcs/BLAKE3/blake3_avx512_x86-64_unix.S.o.requires

.PHONY : ligero/CMakeFiles/test_merkleTree.dir/requires

ligero/CMakeFiles/test_merkleTree.dir/clean:
	cd "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" && $(CMAKE_COMMAND) -P CMakeFiles/test_merkleTree.dir/cmake_clean.cmake
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/clean

ligero/CMakeFiles/test_merkleTree.dir/depend:
	cd "/home/leeweihan/Desktop/ligero++ -1209/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/leeweihan/Desktop/ligero++ -1209" "/home/leeweihan/Desktop/ligero++ -1209/ligero" "/home/leeweihan/Desktop/ligero++ -1209/build" "/home/leeweihan/Desktop/ligero++ -1209/build/ligero" "/home/leeweihan/Desktop/ligero++ -1209/build/ligero/CMakeFiles/test_merkleTree.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : ligero/CMakeFiles/test_merkleTree.dir/depend

