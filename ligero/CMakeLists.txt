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

add_executable(test_zk-1 tests/test_zk-1.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_zk-1 ligero)
target_link_libraries(test_zk-1  -lpthread -lm)

add_executable(test_zk-2 tests/test_zk-2.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_zk-2 ligero)
target_link_libraries(test_zk-2  -lpthread -lm)

add_executable(test_zk-3 tests/test_zk-3.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_zk-3 ligero)
target_link_libraries(test_zk-3  -lpthread -lm)

add_executable(test_zk-4 tests/test_zk-4.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_zk-4 ligero)
target_link_libraries(test_zk-4  -lpthread -lm)

add_executable(test_ipa-1 tests/test_ipa-1.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_ipa-1 ligero)
target_link_libraries(test_ipa-1  -lpthread -lm)

add_executable(test_ipa-2 tests/test_ipa-2.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_ipa-2 ligero)
target_link_libraries(test_ipa-2  -lpthread -lm)

add_executable(test_ipa-3 tests/test_ipa-3.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_ipa-3 ligero)
target_link_libraries(test_ipa-3  -lpthread -lm)

add_executable(test_ipa-4 tests/test_ipa-4.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_ipa-4 ligero)
target_link_libraries(test_ipa-4  -lpthread -lm)

add_executable(test_linear-1 tests/test_linear-1.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_linear-1 ligero)
target_link_libraries(test_linear-1  -lpthread -lm)

add_executable(test_linear-2-1 tests/test_linear-2-1.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_linear-2-1 ligero)
target_link_libraries(test_linear-2-1  -lpthread -lm)

add_executable(test_linear-2-2 tests/test_linear-2-2.cpp
        bcs/BLAKE3/blake3.c
        bcs/BLAKE3/blake3_dispatch.c
        bcs/BLAKE3/blake3_portable.c
        bcs/BLAKE3/blake3_sse2_x86-64_unix.S
        bcs/BLAKE3/blake3_sse41_x86-64_unix.S
        bcs/BLAKE3/blake3_avx2_x86-64_unix.S
        bcs/BLAKE3/blake3_avx512_x86-64_unix.S)
target_link_libraries(test_linear-2-2 ligero)
target_link_libraries(test_linear-2-2  -lpthread -lm)

ENABLE_LANGUAGE(ASM)