# ligerolight
an optimized zero-knowledge argument for blockchain scalability

Ligerolight is a zero-knowledge **S**calable **T**ransparent **AR**gument of **K**nowledge (zk-STARK). It is plausible post-quantum secure, transparent and has poly-lograthimic communication complexity.

**WARNING**: This is an academic proof-of-concept prototype, and in particular has not received careful code review.
This implementation is NOT ready for production use.

=======

The implementation includes:

- An efficient FRI implementation. See ligero/protocols/ldt/fri/fri_ldt.hpp for details and ligero/tests/test_fri.cpp for tests. We refer to some auxiliary functions from [libiop](https://github.com/scipr-lab/libiop/tree/master/libiop/protocols/ldt/fri/fri_aux.hpp). We use the conjectured soundness in [[BCI+20]](https://eprint.iacr.org/2020/654), which is used in [estark]() and [plonky2](https://github.com/mir-protocol/plonky2/blob/main/plonky2/plonky2.pdf). The localization parameter is adjustable.

- An efficient field F_p where p = 2^64 - 2^32 + 1. The FRI protocol requires that there exists multiplicative cosets with order 2^k for large enough k. Our field has a multiplicative coset with size 2^32. Compared with galois field and other prime field in libiop, our field are more efficient to support faster FFTs and multiplication. See the paper for details. The implemented field is in depends/libff/libff/algebra/fields/prime_base/fp_64.tcc.

- A batch zero-knowledge inner product argument. The batch zk-IPA allows proving multiple inner product relations at one time. The communication complexity of this batch zk-IPA is poly-logarithmic to the vector length and the verifier complexity is also logarithmic if using an interactive proof called GKR protocol as delegation. The delegation and FFT circuit refer to [Virgo](https://github.com/sunblaze-ucb/Virgo). The implemented zk-IPA is in ligero/tests/test_PCS.cpp

- Ligerolight. See ligero/tests/test_debugg.cpp. The statement is 
Merkle tree with SHA-256 functions. The circuit is written using compilers in [libsnark](https://github.com/scipr-lab/libsnark). See ligero/gadgetlib1 for detials.
The RS code rate for row and column, the localization parameter array, the circuit size, the security parameter are all adjustable.
We use blake3 hash function for merkle tree, where other compiler languages such as assembly and rust are involved. The test for merkle tree is in ligerolight/ligero/tests/test_merkleTree.cpp.

## Usage

To run the code, first install depends

```bash
sudo apt-get install build-essential cmake git libgmp3-dev libprocps4-dev libboost-all-dev libssl-dev libsodium-dev --fix-missing
git submodule init && git submodule update
```

then build
```bash
mkdir build
cd build
cmake ..
make
```
