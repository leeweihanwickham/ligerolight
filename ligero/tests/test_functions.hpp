//
// Created by wx on 22-10-25.
//

#ifndef LIBIOP_TEST_FUNCTIONS_HPP
#define LIBIOP_TEST_FUNCTIONS_HPP

/* test aurora part */
#include <cstdint>
#include <stdexcept>

#include "ligero/gadgetlib1/gadget.hpp"
#include "ligero/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp"
#include "ligero/gadgetlib1/gadgets/merkle_tree/merkle_tree_check_read_gadget.hpp"



namespace ligero{

/* generate a simple sha256 constraints for */
    template<typename FieldT, typename HashT>
    void generate_sha256_example(protoboard<FieldT> &pb, bool pad = true)
    {
        pb_variable_array<FieldT> primary_pad;
        pb_variable_array<FieldT> variable_pad;

        /* 按照顺序 */
        digest_variable<FieldT> output(pb, SHA256_digest_size, "output");
        primary_pad.allocate(pb, 255, "primary_pad");

        digest_variable<FieldT> left(pb, SHA256_digest_size, "left");
        digest_variable<FieldT> right(pb, SHA256_digest_size, "right");

        HashT sha256(pb,left,right,output,"sha256");
        sha256.generate_r1cs_constraints();
        pb.set_input_sizes(511); // num_inputs = 2^n -1
        const std::size_t constraint_dim = libff::log2(pb.num_constraints());
        const size_t num_constraints = 1ul << constraint_dim;
        const size_t num_variables = (1ul << (libff::log2(pb.num_variables()))) -1;
        variable_pad.allocate(pb, num_variables - pb.num_variables(), "variable_pad");

        const libff::bit_vector hash_bv = libff::int_list_to_bits({0xeffd0b7f, 0x1ccba116, 0x2ee816f7, 0x31c62b48, 0x59305141, 0x990e5c0a, 0xce40d33d, 0x0b1167d1}, 32);
        const libff::bit_vector left_bv = libff::int_list_to_bits({0x426bc2d8, 0x4dc86782, 0x81e8957a, 0x409ec148, 0xe6cffbe8, 0xafe6ba4f, 0x9c6f1978, 0xdd7af7e9}, 32);
        const libff::bit_vector right_bv = libff::int_list_to_bits({0x038cce42, 0xabd366b8, 0x3ede7e00, 0x9130de53, 0x72cdf73d, 0xee825114, 0x8cb48d1b, 0x9af68ad0}, 32);
        left.generate_r1cs_witness(left_bv);
        right.generate_r1cs_witness(right_bv);
        sha256.generate_r1cs_witness();

        pb.set_constraint_sizes(num_constraints);

    }

    template<typename FieldT, typename HashT>
    void generate_merkle_check_read_example(protoboard<FieldT> &pb, const std::size_t tree_depth)
    {
        /* prepare test */
        const size_t digest_len = HashT::get_digest_len();
        std::vector<merkle_authentication_node> path(tree_depth);

        libff::bit_vector prev_hash(digest_len);
        std::generate(prev_hash.begin(), prev_hash.end(), [&]() { return std::rand() % 2; });
        libff::bit_vector leaf = prev_hash;

        libff::bit_vector address_bits;

        /* generate random address and hash values.*/
        size_t address = 0;
        for (long level = tree_depth-1; level >= 0; --level)
        {
            const bool computed_is_right = (std::rand() % 2);
            address |= (computed_is_right ? 1ul << (tree_depth-1-level) : 0);
            address_bits.push_back(computed_is_right);
            libff::bit_vector other(digest_len);
            std::generate(other.begin(), other.end(), [&]() { return std::rand() % 2; });

            libff::bit_vector block = prev_hash;
            block.insert(computed_is_right ? block.begin() : block.end(), other.begin(), other.end()); //拼接计算
            libff::bit_vector h = HashT::get_hash(block);

            path[level] = other; //存储着每层其他的数值

            prev_hash = h;
        }
        libff::bit_vector root = prev_hash;

        /* execute test */
        digest_variable<FieldT> root_digest(pb, digest_len, "output_digest");
        /* pad input*/
        std::size_t target_input_num = (digest_len << 1) - 1; // assume that all the digest_len = 2^k
        pb_variable_array<FieldT> input_pad;
        input_pad.allocate(pb, target_input_num - digest_len, "input_pad");


        digest_variable<FieldT> leaf_digest(pb, digest_len, "input_block");
        pb_variable_array<FieldT> address_bits_va;
        address_bits_va.allocate(pb, tree_depth, "address_bits");
        merkle_authentication_path_variable<FieldT, HashT> path_var(pb, tree_depth, "path_var");
        merkle_tree_check_read_gadget<FieldT, HashT> ml(pb, tree_depth, address_bits_va, leaf_digest, root_digest, path_var, ONE, "ml");

        pb.set_input_sizes(target_input_num);
        /* pad variable*/
        //std::size_t target_variable_num = (1ul << libff::log2(pb.num_variables())) -1;
        //std::size_t target_variable_num = (pb.num_variables()) -1;
        pb_variable_array<FieldT> variable_pad;
        // nei cun wen ti, 2022-11.07
        //variable_pad.allocate(pb, target_variable_num - pb.num_variables(), "variable_pad");

        /* generate  & pad constraint */
        path_var.generate_r1cs_constraints();
        ml.generate_r1cs_constraints();
        //const size_t num_constraints = 1ul << (libff::log2(pb.num_constraints()));
        const size_t num_constraints = pb.num_constraints();
        pb.set_constraint_sizes(num_constraints);

        address_bits_va.fill_with_bits(pb, address_bits);
        assert(address_bits_va.get_field_element_from_bits(pb).as_ulong() == address);
        leaf_digest.generate_r1cs_witness(leaf);
        path_var.generate_r1cs_witness(address, path);
        ml.generate_r1cs_witness();

        /* make sure that read checker didn't accidentally overwrite anything */
        address_bits_va.fill_with_bits(pb, address_bits);
        leaf_digest.generate_r1cs_witness(leaf);
        root_digest.generate_r1cs_witness(root);
        assert(pb.is_satisfied());
    }

}


#endif //LIGERO_TEST_FUNCTIONS_HPP