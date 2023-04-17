add_library(
  ligero

  common/common.cpp
  
  bcs/hashing/blake2b.cpp
  protocols/ldt/ldt_reducer.cpp
  iop/utilities/batching.cpp
  algebra/utils.cpp
)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(sodium)

target_link_libraries(
  ligero

  sodium

  ff
)

target_include_directories(
 ligero

 PUBLIC
 ../depends/libfqfft
 ../depends/libff
)

add_executable(test_debugg tests/test_debugg.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_debugg ligero)
target_link_libraries(test_debugg  -lpthread -lm)

add_executable(test_PCS tests/test_PCS.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_PCS ligero)
target_link_libraries(test_PCS  -lpthread -lm)

add_executable(test_PCS_trivial
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_PCS_trivial ligero)
target_link_libraries(test_PCS_trivial  -lpthread -lm)

add_executable(test_merkleTree tests/test_merkleTree.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_merkleTree ligero)
target_link_libraries(test_merkleTree  -lpthread -lm)

add_executable(test_fri tests/test_fri.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_fri ligero gtest_main)

ENABLE_LANGUAGE(ASM)