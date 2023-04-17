#include <cstdint>
#include <stdexcept>
#include <ctime>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/fields/prime_base/fields_64.hpp>
#include "ligero/algebra/fft.hpp"
#include "ligero/algebra/polynomials/polynomial.hpp"
#include "ligero/protocols/ligero_iop.hpp"
#include <sys/time.h>

#include "ligero/tests/test_functions.hpp"

using namespace ligero;
//std::clock_t start,end,total_start,total_end;

// TODO: reduce a blinding matrix for w
// TODO: reduce the oracle matrix height related to localication_arrary, This is related to both time and proof size
// TODO:

int main(){
    struct timeval start,end,total_start;
//    libff::alt_bn128_pp::init_public_params();
//    typedef libff::alt_bn128_Fr FieldT;
    typedef libff::Fields_64 FieldT;
    const std::size_t RS_extra_dimensions = 2;//2 or 3
    const std::size_t RS_col_extra_dimensions = 3;
    const std::vector<std::size_t> localization_parameter_array({1,4,4});
    const std::size_t security_parameter = 100;//P(soundness error)<2^(-security_parameter)
    const bool make_zk = true;
    const float height_width_ratio = 1ull << 14;
    const std::size_t parallel_opt=0;// 并行选项 0则不并行 其他为并行线程数

    gettimeofday(&total_start, nullptr);
    typedef sha256_two_to_one_hash_gadget<FieldT> HashT;
    protoboard<FieldT> pb;
    // only one sha-256
    //generate_sha256_example<FieldT,HashT>(pb);
    // #hash = tree_depth
    generate_merkle_check_read_example<FieldT,HashT>(pb,8);

    r1cs_constraint_system<FieldT> cs = pb.get_constraint_system();
    //EXPECT_TRUE(pb.is_satisfied());

    gettimeofday(&end, nullptr);
    const float circuit_generate = (end.tv_usec-total_start.tv_usec)/1000000.0 + end.tv_sec-total_start.tv_sec;

    std::map<std::size_t,FieldT> public_input;//z
    public_input.insert(std::pair<std::size_t, FieldT>(pb.primary_input().size()+pb.auxiliary_input().size(), pb.auxiliary_input().back()));

    iop_protocol<FieldT> IOP;
    std::size_t num_variables = cs.num_variables();
    std::size_t num_constraints = cs.num_constraints();

    ligero_iop_parameters<FieldT> params(security_parameter,
                                                 RS_extra_dimensions,
                                                 RS_col_extra_dimensions,
                                                 localization_parameter_array,
                                                 height_width_ratio,
                                                 make_zk,
                                                 multiplicative_coset_type,
                                                 num_constraints,
                                                 num_variables,
                                                 parallel_opt);

    ligero_iop<FieldT> proto(IOP,
                             cs,
                             params);

    gettimeofday(&end, nullptr);
    const float preprocessor = (end.tv_usec-total_start.tv_usec)/1000000.0 + end.tv_sec-total_start.tv_sec;

    gettimeofday(&start, nullptr);
    proto.register_interactions();
    IOP.seal_interaction_registrations();
    proto.register_queries();
    IOP.seal_query_registrations();

    const float mt_arrange_time = proto.produce_oracle(pb.primary_input(), pb.auxiliary_input(), public_input);
    gettimeofday(&end, nullptr);
    const float oracle_generate = (end.tv_usec-start.tv_usec)/1000000.0 + end.tv_sec-start.tv_sec - mt_arrange_time;

    gettimeofday(&start, nullptr);
    const float random_chuanzhi_time=proto.produce_proof(public_input);
    gettimeofday(&end, nullptr);
    const float proof_generate = (end.tv_usec-start.tv_usec)/1000000.0 + end.tv_sec-start.tv_sec-random_chuanzhi_time;

    gettimeofday(&start, nullptr);
    bool correctness = proto.verifier_predicate(public_input);
    gettimeofday(&end, nullptr);
//
    const float verify_time = (end.tv_usec-start.tv_usec)/1000000.0 + end.tv_sec-start.tv_sec;
    const float total_time = (end.tv_usec-total_start.tv_usec)/1000000.0 + end.tv_sec-total_start.tv_sec;
    if(correctness){
        std::cout<<"protocol run correctly!"<<std::endl;
        std::cout<<"circuit_generate_time"<<std::endl;
        std::cout<<circuit_generate<<" s"<<std::endl;
        std::cout<<"preprocessor_time"<<std::endl;
        std::cout<<preprocessor<<" s"<<std::endl;
        std::cout<<"oracle_generate_time"<<std::endl;
        std::cout<<oracle_generate<<" s"<<std::endl;
        std::cout<<"proof_generate_time"<<std::endl;
        std::cout<<proof_generate<<" s"<<std::endl;
        std::cout<<"verify_time"<<std::endl;
        std::cout<<verify_time<<"s"<<std::endl;
        std::cout<<"total_time"<<std::endl;
        std::cout<<total_time<<"s"<<std::endl;
    }
    else{
        std::cout<<"error occurs!"<<std::endl;
    }

    return 0;

}

