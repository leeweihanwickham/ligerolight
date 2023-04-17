//
// Created by wx on 22-10-20.
//

#ifndef LIBIOP_MERKLE_TREE_CHECK_HPP
#define LIBIOP_MERKLE_TREE_CHECK_HPP

#include <ligero/common/data_structures/merkle_tree.hpp>
#include <ligero/gadgetlib1/gadget.hpp>
#include <ligero/gadgetlib1/gadgets/hashes/hash_io.hpp>

namespace ligero{

template<typename FieldT, typename HashT>
class merkle_circuit: public gadget<FieldT>{
private:

    std::vector<HashT> hashers;
    // std::vector<block_variable<FieldT> > hasher_inputs;
    std::vector<digest_variable<FieldT> > internal_output;

    // std::shared_ptr<digest_variable<FieldT> > computed_root; 合并进internal_output里了,是internal_output[0]
    // std::shared_ptr<bit_vector_copy_gadget<FieldT> > check_root;

public:

    const size_t digest_size;
    const size_t tree_depth;
    const size_t leaves_number;
    // 根据实际hash函数的digest_size确定后
    digest_variable<FieldT> root_digest;
    std::vector<digest_variable<FieldT>> leaves_digest;


    merkle_circuit(protoboard<FieldT> &pb,
                    const size_t tree_depth,
                    const digest_variable<FieldT> &root_digest,
                    const std::vector<digest_variable<FieldT>> &leaves_digest, 
                    const std::string &annotation_prefix="");

    void generate_r1cs_constraint();
    void generate_r1cs_witness();

    static libff::bit_vector get_root(const size_t tree_depth,const std::vector<libff::bit_vector> &leaves_bits);
};

template<typename FieldT, typename HashT>
void test_merkle_tree_gadget();

}

#include <ligero/gadgetlib1/gadgets/merkle_tree/merkle_tree_check.tcc>
#endif //LIBIOP_MERKLE_TREE_CHECK_HPP
