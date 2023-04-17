#include <algorithm>
#include <cstdint>
#include <iostream>

#include <gtest/gtest.h>

#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/fields/prime_base/fields_64.hpp>
#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/field_utils/field_utils.hpp>
#include "ligero/algebra/fft.hpp"
#include "ligero/algebra/field_subset/subgroup.hpp"
#include <libff/common/utils.hpp>
#include "ligero/iop/iop.hpp"
#include "ligero/protocols/ldt/fri/fri_ldt.hpp"
#include "ligero/algebra/polynomials/polynomial.hpp"

/*commitment related*/

#include "ligero/bcs/hashing/blake2b.hpp"
#include "ligero/bcs/hashing/dummy_algebraic_hash.hpp"
#include "ligero/bcs/hashing/hashing.hpp"
#include "ligero/bcs/merkle_tree.hpp"
#include "ligero/bcs/Newmerkle.hpp"

namespace ligero {


template<typename FieldT>
bool run_test(const std::size_t codeword_domain_dim,
              const std::vector<std::size_t> localization_parameter_array,
              const std::size_t RS_extra_dimensions,
              const std::size_t poly_degree_bound,
              const std::size_t num_interactive_repetitions = 0, /* auto filled for 129 bits of security if 0 */
              const std::size_t num_query_repetitions = 0)
{

    /**NOTE: must be domain (1<<size, FieldT(1<<size))
     * do not know why now...**/
    field_subset<FieldT> domain(1 << codeword_domain_dim, FieldT(1 << codeword_domain_dim));
    std::size_t codeword_domain_size = 1 << codeword_domain_dim;
    const polynomial<FieldT> poly = polynomial<FieldT>::random_polynomial(poly_degree_bound);
    // a null value
    FRI_verifier<FieldT> *verifier = new FRI_verifier<FieldT>(poly_degree_bound, localization_parameter_array, domain);
    // 3 kinds function, compute the first inteplot
    FRI_prover<FieldT> *prover = new FRI_prover<FieldT>(poly, localization_parameter_array, verifier, domain);
    //std::vector<std::size_t> query_set = {19260817, 20221130, 238456983, 7056978, 12312, 54976};
    std::vector<std::size_t> query_set;
    for (std::size_t i = 0; i < 10; ++i)
    {
        bool is_repeat = true;
        std::size_t val;
        while(is_repeat){
            //val = std::rand() % (codeword_domain_size);
            val = std::rand() % (codeword_domain_size >> localization_parameter_array[0]);
            std::vector<std::size_t>::iterator it;
            it = find(query_set.begin(),query_set.end(),val);
            if(it == query_set.end()){
                is_repeat = false;
            }
        }
        query_set.emplace_back(val);
    }
    std::sort(query_set.begin(),query_set.end());
    for (std::size_t i = 0; i < query_set.size(); i ++)
    {
        query_set[i] += (codeword_domain_size - 1);
    }
    prover->prove(query_set);
    prover->query(query_set);
    return verifier->verify(query_set, prover, prover->final_poly_coeffs);
}

template<typename FieldT>
bool test_inner_product(const field_subset<FieldT> &compute_domain,
                        std::vector<polynomial<FieldT>> IPA_pub_polys,
                        std::vector<polynomial<FieldT>> IPA_sec_polys,
                        FieldT target_sum,
                        std::size_t security_parameter,
                        std::size_t IPA_RS_extra_dimensions,
        // for FRI, need to padding;
        // the max degree for degree{IPA_pub_polys[i]*IPA_sec_polys[i]}_{i}
        // 4H - 3
                        std::size_t poly_degree_bound,
        // the max poly degree for testing polynomial to sumcheck
        // next_power_of_two(poly_degree_bound)
        // 4 H
                        std::size_t summation_degree_bound,
        // a parameter for FRI
                        std::vector<std::size_t> localization_parameter_array) {

    // the summation domain, inner product domain
    const std::size_t compute_domain_dim = libff::log2(compute_domain.num_elements());
    // after encoding
    const std::size_t codeword_domain_dim = compute_domain_dim + IPA_RS_extra_dimensions;
    /** NOTE: must be codeword_domain (1<<size, FieldT(1<<size))
     * do not know why now...**/
    // vertical codeword domain
    field_subset<FieldT> ldt_domain(1 << codeword_domain_dim, FieldT(1 << codeword_domain_dim));
    const std::size_t codeword_domain_size = 1 << codeword_domain_dim;

    // compute parameters
    const long double field_size_bits = (long double)(libff::soundness_log_of_field_size_helper<FieldT>(FieldT::zero()));

    /** [2^{(field_size_bits)} ^ {-e}] * poly_degree_bound = 2^{- security_parameter+1}
     * [2^{(field_size_bits)} ^ {-e}]   = 2^{- security_parameter+1} / poly_degree_bound
     * -e * [(field_size_bits)    = {- security_parameter+1} -  log_2(poly_degree_bound)
     * e      = ({security_parameter-1} +  log_2(poly_degree_bound))/ (field_size_bits)**/
    std::size_t inter_repetition = ceil(( security_parameter - 1 + libff::log2(summation_degree_bound) ) / field_size_bits);
    const std::size_t inter_repetition_num = std::max<size_t> (1, inter_repetition);

    /** IPA_RS_extra_dimensions * query_repetition_num > security_parameter **/
    const std::size_t query_repetition_num = ceil((security_parameter/IPA_RS_extra_dimensions)) + 1;

    std::vector<std::vector<FieldT>> IPA_sec_evaluations;
    IPA_sec_evaluations.resize(IPA_pub_polys.size());
    for (std::size_t i = 0; i < IPA_pub_polys.size(); i ++)
    {
        IPA_sec_evaluations[i] = FFT_over_field_subset(IPA_sec_polys[i].coefficients(), ldt_domain);
    }

    // generate query_Set
    std::vector<std::size_t> IPA_query_set;
    IPA_query_set.clear();
    for (std::size_t i = 0; i < query_repetition_num; ++i) {
        bool is_repeat = true;
        std::size_t val;
        while (is_repeat) {
            val = std::rand() % (codeword_domain_size >> localization_parameter_array[0]);
            std::vector<std::size_t>::iterator it;
            it = find(IPA_query_set.begin(), IPA_query_set.end(), val);
            if (it == IPA_query_set.end()) {
                is_repeat = false;
            }
        }
        IPA_query_set.emplace_back(val);
    }

    for (auto &i: IPA_query_set) {
        i %= (codeword_domain_size >> localization_parameter_array[0]);
    }
    std::vector<std::size_t> tmp;
    std::sort(IPA_query_set.begin(), IPA_query_set.end());
    tmp.push_back(IPA_query_set[0]);
    for (std::size_t i = 1; i < IPA_query_set.size(); i++) {
        if (IPA_query_set[i] != IPA_query_set[i - 1]) {
            tmp.push_back(IPA_query_set[i]);
        }
    }
    IPA_query_set = tmp;

    for (auto &i: IPA_query_set) {
        i += ((codeword_domain_size >> localization_parameter_array[0]) - 1);
    }

    /** It only needs to commit every secret evaluations once as the verifier can construct virtual oracles
     * There is an optimization that every evaluation can be spilt and aggregated together,
     * this is related to the localization array**/

    std::vector<std::vector<FieldT>> commit_matrix;
    commit_matrix.resize((1ull << localization_parameter_array[0]) * (IPA_pub_polys.size() + 1));
    for (std::size_t i = 0; i < IPA_pub_polys.size(); i++) {
        for (std::size_t j = 0; j < (1ull << localization_parameter_array[0]); j++) {
            for (std::size_t k = 0; k < (codeword_domain_size >> localization_parameter_array[0]); k++) {
                commit_matrix[(i * (1ull << localization_parameter_array[0])) + j].push_back
                        (IPA_sec_evaluations[i][k + j * (codeword_domain_size
                                >> localization_parameter_array[0])]);
            }
        }
    }

    std::shared_ptr<ligero::merkle<FieldT>> secret_vector_tree;
    secret_vector_tree.reset(new merkle<FieldT>(
            codeword_domain_size >> localization_parameter_array[0],
            IPA_query_set,
            true
    ));

    ligero::merkleTreeParameter par_for_secret_vector;
    par_for_secret_vector = secret_vector_tree->create_merklePar_of_matrix(commit_matrix);


    // The FFT evluation for secret poly
    std::vector<std::vector<FieldT>> IPA_sec_polys_evluation_on_codeword_domain = IPA_sec_evaluations;

    std::vector<std::vector<FieldT>> IPA_pub_polys_evluation_on_codeword_domain;
    IPA_pub_polys_evluation_on_codeword_domain.resize(IPA_pub_polys.size());
    for (std::size_t i = 0; i < IPA_sec_polys.size() ; i++) {
        IPA_pub_polys_evluation_on_codeword_domain[i] = (FFT_over_field_subset(IPA_pub_polys[i].coefficients(), ldt_domain));
    }

    std::size_t padding_degree = poly_degree_bound - summation_degree_bound;

    /** TODO: the input of verifier can be IPA_pub_polys_evaluation_on_codeword_domain
     * But now is poly itself
     * The actual use is evaluation_at_point
     * Need comparison**/
    libff::enter_block("Setting Inner Product Verifier");
    Inner_product_verifier<FieldT> verifier(std::move(IPA_pub_polys), compute_domain,
                                            std::move(IPA_sec_polys_evluation_on_codeword_domain),
                                            std::move(IPA_pub_polys_evluation_on_codeword_domain),
                                            padding_degree, poly_degree_bound, localization_parameter_array,
                                            ldt_domain, target_sum, inter_repetition_num);
    libff::leave_block("Setting Inner Product Verifier");

    libff::enter_block("Inner Product Prover");
    libff::enter_block("Setting Inner Product Prover and compute the first round");
    // TODO: There is no need to commit v_trees in the prover
    // TODO: The commitment for evaluation can also be split
    Inner_product_prover<FieldT> prover(std::move(IPA_pub_polys), std::move(IPA_sec_polys), IPA_query_set,
                                        localization_parameter_array, poly_degree_bound, verifier, ldt_domain,
                                        inter_repetition_num);
    libff::leave_block("Setting Inner Product Prover and compute the first round");
    libff::enter_block("Proving all the remained rounds for FRI");
    prover.prove(IPA_query_set);
    libff::leave_block("Proving all the remained rounds for FRI");
    libff::leave_block("Inner Product Prover");

    /** Note that the most time for verifier is in verify the first round
     * where Merkle hash tree accouts about half
     * construct f(q) accouts about half **/
    libff::enter_block("Inner Product Verifier");
    bool result = verifier.verify(IPA_query_set, &prover);
    bool result1 = secret_vector_tree->verify_merkle_commit(par_for_secret_vector);
    if (result && result1){
        libff::print_indent(); printf("Protocol runs successfully! \n");
    }
    else {
        libff::print_indent(); printf("error occurs! \n");
    }
    libff::leave_block("Inner Product Verifier");


    /**Proof size compute**/
    /** NOTE: for v_trees_hashes and h_tree_hashes
     * search TODO: proof size: v_trees related
     * and TODO: proof size: h_tree related
     * and TODO: proof size: FRI_trees related **/
    std::size_t v_trees_hashes = par_for_secret_vector.path_lenth;
    std::size_t h_tree_hashes = prover.h_tree_lenth;
    std::size_t FRI_trees_hashes = prover.FRI_tree_lenth ;

    /** in Ligerolight, there will only be 1 tree for sec_polys and v_trees_hashes **/
    //TODO the root number needs to be changed as the column width
    std::size_t proof_size_hash_number = IPA_sec_polys.size() + 1 + (localization_parameter_array.size() - 2) // Merkle tree roots number
                                         + v_trees_hashes  + h_tree_hashes + FRI_trees_hashes; // Merkle authenitacation path number.

    std::size_t hash_size = 256;
    std::size_t proof_size_hash = (proof_size_hash_number * hash_size)/1024/8;

    std::size_t proof_size_field_number = 0;
    std::size_t degree_decrease = 0 ;
    for (std::size_t i = 0 ; i < localization_parameter_array.size() ; i++)
    {
        if (i==0)
        {
            proof_size_field_number +=  inter_repetition_num * IPA_sec_polys.size() * ( 1 << localization_parameter_array[i]);
            degree_decrease += ( 1 << localization_parameter_array[i]);
        }
        else{
            proof_size_field_number +=  inter_repetition_num * ( 1 << localization_parameter_array[i]);
            degree_decrease += ( 1 << localization_parameter_array[i]);
        }
    }
    // The final_poly proof size
    proof_size_field_number += inter_repetition_num * (poly_degree_bound - degree_decrease + 1) ;
    std::size_t proof_size_field = (proof_size_field_number * field_size_bits)/1024/8;

    const std::size_t proof_size = proof_size_field + proof_size_hash;
    std::cout << "proof_size_field is " << proof_size_field << std::endl;
    std::cout << "proof_size_hash is " << proof_size_hash << std::endl;

    return result;
}


TEST(FRIMultiplicativeTrueTest, SimpleTest) {
//    libff::alt_bn128_pp::init_public_params();
//
//    typedef libff::alt_bn128_Fr FieldT;

//    libff::edwards_pp::init_public_params();
//    typedef libff::edwards_Fr FieldT;

    typedef libff::Fields_64 FieldT;

    /* Common parameters */
    const std::size_t codeword_domain_dim = 12;
    const std::size_t RS_extra_dimensions = 3; /* \rho = 2^{-RS_extra_dimensions} */

    std::vector<std::size_t> localization_parameter_array(4, 2);
    localization_parameter_array.insert(localization_parameter_array.begin(), 1);

    const std::size_t poly_degree_bound = 1ull << (codeword_domain_dim - RS_extra_dimensions);

    bool result = run_test<FieldT>(codeword_domain_dim, localization_parameter_array, RS_extra_dimensions, poly_degree_bound);
    EXPECT_TRUE(result);
}

TEST(InnerProductTest, SimpleTest) {
//        libff::alt_bn128_pp::init_public_params();
//
//        typedef libff::alt_bn128_Fr FieldT;
    typedef libff::Fields_64 FieldT;

    /* Common parameters */
    const std::size_t codeword_domain_dim = 11;
    const std::size_t RS_extra_dimensions = 2; /* \rho = 2^{-RS_extra_dimensions} */
    const std::size_t instance = 10;

    std::vector<std::size_t> localization_parameter_array({3,2,2});
    //localization_parameter_array.insert(localization_parameter_array.begin(), 1);

    const std::size_t poly_degree_bound = 1ull << (codeword_domain_dim - RS_extra_dimensions);
    field_subset<FieldT> compute_domain(poly_degree_bound>>1, FieldT::random_element());


    std::vector<polynomial<FieldT>> IPA_pub_polys;
    IPA_pub_polys.resize(instance);
    for (std::size_t i = 0; i < IPA_pub_polys.size() ; i++)
    {
        IPA_pub_polys[i] = polynomial<FieldT>::random_polynomial((poly_degree_bound >>1 ));
//        std::vector<FieldT> coefficient = IPA_pub_polys[i].coefficients();
//        while (true) {
//            if (coefficient.back() != FieldT(0)) {
//                break;
//            }
//            coefficient.pop_back();
//        }
    }

    std::vector<polynomial<FieldT>> IPA_sec_polys;
    IPA_sec_polys.resize(instance);
    for (std::size_t i = 0; i < IPA_sec_polys.size() ; i++)
    {
        IPA_sec_polys[i] = polynomial<FieldT>::random_polynomial((poly_degree_bound >>1) +1);

    }


    FieldT target_sum = FieldT::zero();
    // add
    for (auto &i: compute_domain.all_elements()) {
        for (std::size_t j = 0; j < IPA_pub_polys.size(); j++ )
        {
            target_sum += IPA_pub_polys[j].evaluation_at_point(i) * IPA_sec_polys[j].evaluation_at_point(i);
        }
    }

    // compute_domain : systematic domain vertical
    bool result = test_inner_product<FieldT>(compute_domain, IPA_pub_polys, IPA_sec_polys, target_sum, 100, RS_extra_dimensions,
                                             poly_degree_bound, poly_degree_bound, localization_parameter_array);
    EXPECT_TRUE(result);
}


}
