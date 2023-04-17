//
// Created by wx on 22-10-20.
//

#ifndef LIBIOP_MERKLE_TREE_CHECK_TCC
#define LIBIOP_MERKLE_TREE_CHECK_TCC

namespace ligero{

template<typename FieldT,typename HashT>
merkle_circuit<FieldT,HashT>::merkle_circuit(protoboard<FieldT> &pb,
                                    const size_t tree_depth,
                                    const digest_variable<FieldT> &root_digest,
                                    const std::vector<digest_variable<FieldT>> &leaves_digest, 
                                    const std::string &annotation_prefix):
    gadget<FieldT>(pb, annotation_prefix),
    tree_depth(tree_depth),
    root_digest(root_digest),
    leaves_digest(leaves_digest),
    digest_size(HashT::get_digest_len()),
    leaves_number(leaves_digest.size())
{
    assert(tree_depth > 0);
    assert(leaves_number == (1ul << tree_depth));
    /* Debug用
    for(size_t i = 0; i < leaves_number; ++i)
    {
        leaves_digest[i] = std::make_shared<digest_variable<FieldT>>(pb,digest_size, 
                                                FMT(this->annotation_prefix, " leaf_%zu", i));
    }
    */

    for(size_t i = 0; i < leaves_number - 1; ++i)
    {
        internal_output.emplace_back(digest_variable<FieldT>(pb, digest_size, FMT(this->annotation_prefix, " internal_output_%zu", i)));
    }
    
    /* 为计算第level层的hash值, level=0表示根所在的层 */
    for(size_t level = 0, node_number = 1; level < tree_depth;  ++level)
    {
        if(level == tree_depth-1 ) {
            std::cout << "叶子" << std::endl;
            for (size_t i = 0; i < node_number; ++i) {
                hashers.emplace_back(HashT(pb, leaves_digest[2 * i], leaves_digest[2 * i + 1],
                                           internal_output[node_number - 1 + i],
                                           FMT(this->annotation_prefix, " load_hashers_%zu", node_number - 1 + i)));
                /*std::cout << "当前绑定的（f,s,s）为（" << node_number - 1 + i << " " << 2 * i  << " " << 2 * i + 1
                          << std::endl;                                             */

            }
        }
        else
        {
            for(size_t i = 0; i < node_number; ++i)
            {
                size_t father = node_number-1+i;
                hashers.emplace_back(HashT(pb, internal_output[2*father+1], internal_output[2*father+2],
                                           internal_output[father],
                                           FMT(this->annotation_prefix, " load_hashers_%zu", node_number-1+i)));
                //std::cout << "当前绑定的（f,s,s）为（" <<  father << " " << 2*father+1 << " " << 2*father+2 << "）" << std::endl;
            }
        }
        
        node_number += node_number;
    }
    
    /*确保我们通过电路计算出的root跟提供的root一致
    这里的read_successful被我用ONE直接代替了，不知道有没有什么后果...*/
    //check_root.reset(new bit_vector_copy_gadget<FieldT>(pb, internal_output[0].bits, root_digest.bits, ONE, FieldT::capacity(), FMT(annotation_prefix, " check_root")));
}

template<typename FieldT,typename HashT>
void merkle_circuit<FieldT,HashT>::generate_r1cs_constraint()
{
    assert(hashers.size() == (leaves_number -1));
    for (size_t i = hashers.size(); i > 0; --i)
    {
        hashers[i-1].generate_r1cs_constraints(false);
    }
    
    //check_root->generate_r1cs_constraints(false, false);
}

template<typename FieldT, typename HashT>
void merkle_circuit<FieldT, HashT>::generate_r1cs_witness()
{
    for (size_t i = hashers.size(); i > 0; --i)
    {
        hashers[i-1].generate_r1cs_witness();
    }

    //check_root->generate_r1cs_witness();
}

template<typename FieldT, typename HashT>
libff::bit_vector merkle_circuit<FieldT,HashT>::get_root(const size_t tree_depth,
                                                         const std::vector<libff::bit_vector> &leaves_bits)
{
    protoboard<FieldT> pb;

    size_t digest_len = HashT::get_digest_len();
    digest_variable<FieldT> root(pb, digest_len, "root");
    std::vector<digest_variable<FieldT>> leaves;
    for(int i = 0; i < leaves_bits.size(); i++)
    {
        leaves.emplace_back(digest_variable<FieldT>(pb, digest_len, FMT("left_digests_%zu", i)));
    }
    merkle_circuit<FieldT,HashT> mc(pb,tree_depth,root,leaves,"mc");
    for(int i= 0; i < leaves_bits.size(); i++)
    {
        leaves[i].generate_r1cs_witness(leaves_bits[i]);
    }
    mc.generate_r1cs_witness();
    return root.get_digest();
}

template<typename FieldT, typename HashT>
void test_merkle_tree_gadget()
{
    const size_t digest_len = HashT::get_digest_len();
    const size_t tree_depth = 4;
    const size_t leaves_number = 1ul << tree_depth;
    // generate leaves
    std::vector<libff::bit_vector> leaves_bits(leaves_number);
    for(int i = 0; i < leaves_number; i++)
    {
        leaves_bits[i].resize(digest_len);
        std::generate(leaves_bits[i].begin(),leaves_bits[i].end(), [&](){ return std::rand()%2; });
    }

    protoboard<FieldT> pb;
    digest_variable<FieldT> root_digest(pb, digest_len, "output_digest");
    std::vector<digest_variable<FieldT>> leaves_digest;
    for(int i = 0; i < leaves_number; i++)
    {
        leaves_digest.emplace_back(digest_variable<FieldT>(pb, digest_len, FMT("left_digests_%zu", i)));
    }
    merkle_circuit<FieldT,HashT> mc(pb,tree_depth,root_digest,leaves_digest,"mc");
    mc.generate_r1cs_constraint();

    for(int i = 0; i < leaves_number; i++)
    {
        leaves_digest[i].generate_r1cs_witness(leaves_bits[i]);
    }

    mc.generate_r1cs_witness();

    assert(pb.is_satisfied());
    /*
    for(int i = 0; i < mc.internal_output.size(); i++)
    {
        libff::bit_vector root_bit = mc.internal_output[i].get_digest();
        for(int j = 0; j < root_bit.size(); j++){
            std::cout << root_bit[j] ;
        }
        std::cout << std::endl;
    }
    */
}

}


#endif //LIBIOP_MERKLE_TREE_CHECK_TCC
