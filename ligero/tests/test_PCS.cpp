//
// Created by leeweihan on 23-4-10.
//
#include <algorithm>
#include <cstdint>
#include <iostream>


#include <libff/algebra/fields/prime_base/fields_64.hpp>
#include <libff/algebra/field_utils/field_utils.hpp>
#include "ligero/algebra/fft.hpp"
#include "ligero/algebra/field_subset/subgroup.hpp"
#include <libff/common/utils.hpp>
#include "ligero/iop/iop.hpp"
#include "ligero/protocols/ldt/fri/fri_ldt.hpp"
#include "ligero/protocols/fft_circuit_GKR.hpp"
#include <sys/time.h>

using namespace ligero;


// We consider the case: (1,x) times (1,x^2)
// According to the order, the latter term size is forever 2
template<typename FieldT>
std::vector<FieldT> TensorProductOfTwo(std::vector<FieldT> vector_1,
                                       std::vector<FieldT> vector_2){
    assert(vector_2.size()==2);
    std::vector<FieldT> result;
    result.resize(vector_1.size() * vector_2.size());
    for (std::size_t i = 0; i < vector_2.size(); i ++)
    {
        for (std::size_t j = 0; j < vector_1.size(); j ++)
        {
            result[i*vector_1.size()+j] = vector_2[i] * vector_1[j];
        }
    }
    return result;
};

template<typename FieldT>
FieldT FieldPow(FieldT x, std::size_t k){
    FieldT result = x;
    for (std::size_t i = 0; i < k; i ++)
    {
        result *= x;
    }
    return result;
};

template<typename FieldT>
std::vector<FieldT> TensorProduct(std::size_t variable_num,
                                  std::size_t degree_dim,
                                  std::vector<FieldT> variable_eva){
    assert(variable_eva.size() == variable_num) ;
    std::vector<FieldT> result;
    std::size_t poly_size = std::size_t (pow ((1ull << degree_dim), variable_num));
    result.resize(poly_size);
    std::vector<FieldT> initial_vector = {1};
    for (std::size_t i = 0; i < variable_num; i ++)
    {
        for (std::size_t j = 0; j < degree_dim; j ++)
        {
            std::vector<FieldT> right_vector = {1, FieldPow(variable_eva[i],j)};
            initial_vector = TensorProductOfTwo(initial_vector,right_vector);
        }
    }
    return initial_vector;
};

int main() {

    typedef libff::Fields_64 FieldT;

//    //test begin
//    std::vector<FieldT> v1 = {1};
//    std::vector<FieldT> v2 = {2,3};
//    std::cout << TensorProduct(2,2,v2);
//    // test end

    const std::size_t repeat_num = 10;

    // common parameters
    const std::size_t variable_num = 10;
    // degree dim and degree bound for every variable
    const std::size_t degree_dim = 1;
    const std::size_t degree = (1ull << degree_dim) - 1;
    // batch number
    const std::size_t instance = 1;
    // rho
    const std::size_t RS_extra_dimension = 3;
    // eta
    std::vector<std::size_t> localization_parameter_array({2,2,2,2});
    // lambda
    const std::size_t security_parameter = 100;
    // |F|
    const std::size_t field_size_bits =
            (long double) (libff::soundness_log_of_field_size_helper<FieldT>(FieldT::zero()));

    libff::enter_block("Configure parameters and generate witness");

    /** Compute other parameters
     * the query repetition parameter according to the conjecture; l
     * the max poly degree when invoking sumcheck instance; k
     * random challenges of verifier;
     * the interactive repetition parameter;
     * the hash output size; |H| TODO: how to set?
     * achieved soundness bits  **/

    const std::size_t poly_size = std::size_t(pow(degree + 1, variable_num));
    std::cout << "poly_size is " << poly_size << std::endl;

    // l - query repetition parameter
    const std::size_t query_repetition_parameter = ceil(double(security_parameter) / RS_extra_dimension);
    // the degree bound of the polynomial to be summed
    const std::size_t sum_degree_bound = 2 * poly_size - 1;
    std::size_t FRI_degree_bound = libff::round_to_next_power_of_2(sum_degree_bound);

    // interactive repetition parameter
    /** the interactive repetition parameter of FRI
     *  [2^{(field_size_bits)} ^ {-e}] * poly_degree_bound = 2^{- security_parameter}
     * [2^{(field_size_bits)} ^ {-e}]   = 2^{- security_parameter} / poly_degree_bound
     * -e * [(field_size_bits)    = {- security_parameter} -  log_2(poly_degree_bound)
     * e      = ({security_parameter} +  log_2(poly_degree_bound))/ (field_size_bits)**/
    std::size_t inter_repetition = ceil(double(security_parameter + libff::log2(FRI_degree_bound)) / field_size_bits);
    const std::size_t inter_repetition_parameter = std::max<size_t>(1, inter_repetition);

    // hash function parameter
    const std::size_t hash_ouput_size = 256;

    /** compute achieved soundness parameter **/
    std::size_t FRI_interactive_error =
            (inter_repetition_parameter * field_size_bits) - ceil(libff::log2(FRI_degree_bound));
    std::size_t FRI_query_error = query_repetition_parameter * RS_extra_dimension;
    std::size_t achieved_soundness = std::min<std::size_t>({FRI_interactive_error, FRI_query_error});


    /** determine the polynomials and domains **/

    field_subset<FieldT> summation_domain(poly_size);
    field_subset<FieldT> codeword_domain(FRI_degree_bound << RS_extra_dimension,
                                         FieldT(FRI_degree_bound << RS_extra_dimension));
    // every instance will have two pairs of polynomials,
    std::size_t poly_number = instance;

    printf("\nPolynomial commitment scheme parameters\n");
    libff::print_indent();
    printf("* target security parameter = %zu\n", security_parameter);
    libff::print_indent();
    printf("* achieved security parameter = %zu\n", achieved_soundness);
    libff::print_indent();
    printf("* variable number = %zu\n", variable_num);
    libff::print_indent();
    printf("* variable degree bound = %zu\n", degree);
    libff::print_indent();
    printf("* RS extra dimensions = %zu\n", RS_extra_dimension);
    libff::print_indent();
    printf("* field size bits = %zu\n", field_size_bits);
    libff::print_indent();
    printf("* FRI interactive repetitions = %zu\n", inter_repetition_parameter);
    libff::print_indent();
    printf("* FRI query repetitions = %zu\n", query_repetition_parameter);
    libff::print_indent();
    printf("* summation degree bound = %zu\n", sum_degree_bound);
    libff::print_indent();
    printf("* FRI degree bound = %zu\n", FRI_degree_bound);

    libff::leave_block("Configure parameters and generate witness");

    libff::enter_block("Generate secret and public vectors");
    /** Configure the public vector
    * Without of generality, we assume every polynomial has its unique evaluation point**/
    std::vector<std::vector<FieldT>> variables_eva;
    variables_eva.resize(instance);

    for (std::size_t i = 0; i < instance; i++) {
        variables_eva[i].resize(variable_num);
        for (std::size_t j = 0; j < variable_num; j++) {
            variables_eva[i][j] = FieldT::random_element();
        }
    }

    std::vector<std::vector<FieldT>> public_vectors;
    public_vectors.resize(instance);
    for (std::size_t i = 0; i < instance; i ++)
    {
        public_vectors[i] = TensorProduct(variable_num, degree_dim, variables_eva[i]);
    }

    /** Generate random secret vectors **/

    std::vector<std::vector<FieldT>> secret_vectors;
    secret_vectors.resize(instance);
    for (std::size_t i = 0; i < instance; i++) {
        secret_vectors[i].resize(poly_size);
        for (std::size_t j = 0; j < poly_size; j++) {
            secret_vectors[i][j] = FieldT::random_element();
        }
    }

    libff::leave_block("Generate secret and public vectors");

    std::size_t repeat;
    struct timeval prover_start, prover_end;
    float prover_time = 0;
    struct timeval verifier_start, verifier_end;
    float verifier_time = 0;
    struct timeval pubpoly_start, pubpoly_end;
    float pubpoly_time = 0;

    double total_proof_size_field = 0;
    double total_proof_size_hash = 0;

    for (repeat = 0; repeat < repeat_num; repeat++) {

        /** PCS Prover **/

        libff::enter_block("PCS prover");
        gettimeofday(&prover_start, nullptr);

        libff::enter_block("Initial polynomials and target sum");

        libff::enter_block("Initial secret polynomials and compute evaluations");

        // constant polynomial
        vanishing_polynomial<FieldT> vanishing_polynomial(summation_domain);

        std::vector<polynomial<FieldT>> IPA_sec_polys;
        IPA_sec_polys.resize(poly_number);

        std::vector<std::vector<FieldT>> IPA_sec_evaluations;
        IPA_sec_evaluations.resize(poly_number);


        for (std::size_t i = 0; i < instance; i++) {
            IPA_sec_polys[i] = polynomial<FieldT>(IFFT_over_field_subset(secret_vectors[i], summation_domain));
            IPA_sec_evaluations[i] = FFT_over_field_subset(IPA_sec_polys[i].coefficients(), codeword_domain);
        }


        libff::leave_block("Initial secret polynomials and compute evaluations");

        libff::enter_block("Initial public polynomials and compute evaluations");

        gettimeofday(&pubpoly_start, nullptr);

        std::vector<polynomial<FieldT>> IPA_pub_polys;
        IPA_pub_polys.resize(instance);

        std::vector<std::vector<FieldT>> IPA_pub_evaluations;
        IPA_pub_evaluations.resize(poly_number);

        std::vector<std::vector<FieldT>> public_poly_evaluations;
        IPA_pub_evaluations.resize(instance);

        for (std::size_t i = 0; i < instance; i++) {
            IPA_pub_polys[i] = polynomial<FieldT>(IFFT_over_field_subset(public_vectors[i], summation_domain));
            IPA_pub_evaluations[i] = FFT_over_field_subset(IPA_pub_polys[i].coefficients(), codeword_domain);
        }


        /** generate \gamma(x), it equals to add a secret poly \gamma(x) and a public poly 1 **/
        polynomial<FieldT> gamma = polynomial<FieldT>::random_polynomial(sum_degree_bound);
        std::vector<FieldT> gamma_eva = FFT_over_field_subset(gamma.coefficients(), codeword_domain);

        std::vector<FieldT> constant_vec(1, FieldT::one());
        polynomial<FieldT> constant_poly = polynomial<FieldT>(std::move(constant_vec));
        std::vector<FieldT> constant_poly_eva = std::vector<FieldT>(codeword_domain.num_elements(), FieldT::one());

        IPA_sec_polys.resize(instance + 1);
        IPA_sec_polys[instance] = gamma;
        IPA_pub_polys.resize(instance + 1);
        IPA_pub_polys[instance] = constant_poly;

        IPA_sec_evaluations.resize(instance + 1);
        IPA_sec_evaluations[instance] = gamma_eva;
        IPA_pub_evaluations.resize(instance + 1);
        IPA_pub_evaluations[instance] = constant_poly_eva;

        gettimeofday(&pubpoly_end, nullptr);

        // compute Gamma
        FieldT Gamma = FieldT::zero();
        std::cout << "IPA_sec_polys.size() is " << IPA_sec_polys.size() << std::endl;
        std::size_t extended_summation_domain_size = libff::round_to_next_power_of_2(sum_degree_bound);
        field_subset<FieldT> extended_summation_domain(extended_summation_domain_size);
        std::vector<FieldT> gamma_eva_on_summation = FFT_over_field_subset(gamma.coefficients(),
                                                                           extended_summation_domain);
        for (std::size_t j = 0; j < summation_domain.num_elements(); j++) {
            std::size_t idx = extended_summation_domain.reindex_by_subset(summation_domain.dimension(), j);
            Gamma += gamma_eva_on_summation[idx] * FieldT::one();
        }

        //std::cout << "Gamma is " << Gamma << std::endl;

        // check the inner product argument
        /**Note here must be next to power of 2**/
        FieldT val1 = FieldT::zero();
        for (std::size_t i = 0 ; i < (IPA_sec_polys.size()); i++){
            std::vector<FieldT> secret_eval = FFT_over_field_subset(IPA_sec_polys[i].coefficients(), extended_summation_domain);
            std::vector<FieldT> public_eval = FFT_over_field_subset(IPA_pub_polys[i].coefficients(), summation_domain);
             //here must be public_eval.size()
            for(std::size_t j = 0 ; j < summation_domain.num_elements() ; j++){
                std::size_t idx = extended_summation_domain.reindex_by_subset(summation_domain.dimension(), j);
                val1+=secret_eval[idx]*public_eval[j];
            }
            // another straightforward way
    //        for (auto j: summation_domain.all_elements())
    //        {
    //            val1 += (IPA_sec_polys[i].evaluation_at_point(j) * IPA_pub_polys[i].evaluation_at_point(j));
    //        }
        }
        FieldT target_sum = val1;

        pubpoly_time += (pubpoly_end.tv_usec-pubpoly_start.tv_usec)/1000000.0 + pubpoly_end.tv_sec-pubpoly_start.tv_sec;
        libff::leave_block("Initial public polynomials and compute evaluations");

        libff::enter_block("GKR prover");

        struct timeval prover_gkr_begin,prover_gkr_end;
        gettimeofday(&prover_gkr_begin, nullptr);

        float p_fft_time = 0;
        float GKR_prove_time = 0;
        size_t fft_GKR_proof_size_field = 0;

        std::vector<FieldT> poly_coe_GKR(IPA_pub_polys[0].coefficients().size(),FieldT::zero());
        std::vector<FieldT> poly_eva_GKR(codeword_domain.num_elements(), FieldT::zero());

        for (std::size_t j = 0; j < IPA_pub_polys[0].coefficients().size(); j ++)
        {
            for (std::size_t i = 0; i < IPA_pub_polys.size(); i ++)
            {
                poly_coe_GKR[j] += IPA_pub_polys[i].coefficients()[j];
            }
        }

        for (std::size_t j = 0; j < IPA_pub_evaluations[0].size(); j ++)
        {
            for (std::size_t i = 0; i < IPA_pub_evaluations.size(); i ++)
            {
                poly_eva_GKR[j] += IPA_pub_evaluations[i][j];
            }
        }

        std::vector<fft_GKR<FieldT>> fft_GKRs;

        fft_GKRs.emplace_back(codeword_domain.coset(),0);

        fft_GKRs[0].prover_compute(poly_coe_GKR,poly_eva_GKR);
        p_fft_time += fft_GKRs[0].p_fft_time;

        std::cout << "p_fft_time is " << p_fft_time << std::endl;
        GKR_prove_time += fft_GKRs[0].p_time;

        std::cout << "GKR_prove_time is " << GKR_prove_time << std::endl;

        fft_GKR_proof_size_field += fft_GKRs[0].proof_size;

        gettimeofday(&prover_gkr_end, nullptr);


        libff::leave_block("GKR prover");

        libff::leave_block("Initial polynomials and target sum");

        libff::enter_block("Generating Merkle tree roots");

        std::vector<std::size_t> IPA_query_set;
        IPA_query_set.clear();
        for (std::size_t i = 0; i < query_repetition_parameter; ++i) {
            bool is_repeat = true;
            std::size_t val;
            while (is_repeat) {
                val = std::rand() % (codeword_domain.num_elements() >> localization_parameter_array[0]);
                std::vector<std::size_t>::iterator it;
                it = find(IPA_query_set.begin(), IPA_query_set.end(), val);
                if (it == IPA_query_set.end()) {
                    is_repeat = false;
                }
            }
            IPA_query_set.emplace_back(val);
        }

        for (auto &i: IPA_query_set) {
            i %= (codeword_domain.num_elements() >> localization_parameter_array[0]);
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
            i += ((codeword_domain.num_elements() >> localization_parameter_array[0]) - 1);
        }

        /** It only needs to commit every secret evaluations once as the verifier can construct virtual oracles
         * There is an optimization that every evaluation can be spilt and aggregated together,
         * this is related to the localization array**/

        std::vector<std::vector<FieldT>> commit_matrix;
        commit_matrix.resize((1ull << localization_parameter_array[0]) * (instance + 1));
        for (std::size_t i = 0; i < instance; i++) {
            for (std::size_t j = 0; j < (1ull << localization_parameter_array[0]); j++) {
                for (std::size_t k = 0; k < (codeword_domain.num_elements() >> localization_parameter_array[0]); k++) {
                    commit_matrix[(i * (1ull << localization_parameter_array[0])) + j].push_back
                            (IPA_sec_evaluations[i][k + j * (codeword_domain.num_elements()
                                    >> localization_parameter_array[0])]);
                }
            }
        }

        for (std::size_t j = 0; j < (1ull << localization_parameter_array[0]); j++) {
            for (std::size_t k = 0; k < (codeword_domain.num_elements() >> localization_parameter_array[0]); k++) {
                commit_matrix[(instance * (1ull << localization_parameter_array[0])) + j].push_back
                        (gamma_eva[k + j * (codeword_domain.num_elements() >> localization_parameter_array[0])]);
            }
        }

        // true is every column put in one leaf
        std::shared_ptr<ligero::merkle<FieldT>> secret_vector_tree;
        secret_vector_tree.reset(new merkle<FieldT>(
                codeword_domain.num_elements() >> localization_parameter_array[0],
                IPA_query_set,
                true
        ));

        ligero::merkleTreeParameter par_for_secret_vector;
        par_for_secret_vector = secret_vector_tree->create_merklePar_of_matrix(commit_matrix);

        libff::leave_block("Generating Merkle tree roots");

        libff::enter_block("Setting parameters");

        // min padding degree
        std::size_t padding_degree = FRI_degree_bound - sum_degree_bound;

        libff::leave_block("Setting parameters");

        libff::enter_block("Setting inner product verifier");

        std::shared_ptr<Inner_product_verifier<FieldT>> IPA_verifier_;
        std::shared_ptr<Inner_product_prover<FieldT>> IPA_prover_;
        std::vector<polynomial<FieldT>> IPA_pub_polys_2 = IPA_pub_polys;

        IPA_verifier_.reset(new Inner_product_verifier<FieldT>(std::move(IPA_pub_polys), summation_domain,
                                                               std::move(IPA_sec_evaluations),
                                                               std::move(IPA_pub_evaluations), padding_degree,
                                                               FRI_degree_bound,
                                                               localization_parameter_array, codeword_domain,
                                                               target_sum, inter_repetition_parameter));

        libff::leave_block("Setting inner product verifier");

        libff::enter_block("Inner Product Prover");
        libff::enter_block("Setting Inner Product Prover and compute the first round");
        IPA_prover_.reset(
                new Inner_product_prover<FieldT>(std::move(IPA_pub_polys_2), std::move(IPA_sec_polys), IPA_query_set,
                                                 localization_parameter_array, FRI_degree_bound, *(IPA_verifier_),
                                                 codeword_domain,
                                                 inter_repetition_parameter));
        libff::leave_block("Setting Inner Product Prover and compute the first round");

        libff::enter_block("Proving all the remained rounds for FRI");
        IPA_prover_->prove(IPA_query_set);
        libff::leave_block("Proving all the remained rounds for FRI");

        libff::leave_block("Inner Product Prover");

        gettimeofday(&prover_end, nullptr);
        prover_time += (prover_end.tv_usec-prover_start.tv_usec)/1000000.0 + prover_end.tv_sec-prover_start.tv_sec - p_fft_time;
        libff::leave_block("PCS prover");

        /** Range Proof Verifier **/

        libff::enter_block("PCS Verifier");
        gettimeofday(&verifier_start, nullptr);

        libff::enter_block("Inner product Verifier");
        bool result1 = IPA_verifier_->verify(IPA_query_set, &(*IPA_prover_));
        libff::leave_block("Inner product Verifier");

        libff::enter_block("Merkle tree Verifier");
        bool result2 = secret_vector_tree->verify_merkle_commit(par_for_secret_vector);
        libff::leave_block("Merkle tree Verifier");


        float verify_fft_GKR_time = 0;
        for(size_t i; i < fft_GKRs.size(); i++)
        {
            bool result = fft_GKRs[i].verifier_predicate();
            verify_fft_GKR_time += fft_GKRs[i].v_time;
        }

        gettimeofday(&verifier_end, nullptr);
        verifier_time += (verifier_end.tv_usec-verifier_start.tv_usec)/1000000.0 + verifier_end.tv_sec-verifier_start.tv_sec;
        libff::leave_block("PCS Verifier");

        if (!(result1 && result2)) {
            libff::print_indent();
            printf("error occurs! \n");
        } else {
            libff::print_indent();
            printf("protocol runs successfully! \n");
        }

        /** Compute the proof size **/

        std::size_t secret_vector_tree_hashes = par_for_secret_vector.path_lenth;
        std::cout << "FRI_tree_hashes 1 " << secret_vector_tree_hashes << std::endl;
        std::size_t h_tree_hashes = IPA_prover_->h_tree_lenth;
        std::cout << "FRI_tree_hashes 2 " << h_tree_hashes << std::endl;
        std::size_t FRI_trees_hashes = IPA_prover_->FRI_tree_lenth;
        std::cout << "FRI_tree_hashes 3 " << FRI_trees_hashes << std::endl;
        std::size_t proof_size_path_hash_number = secret_vector_tree_hashes + h_tree_hashes + FRI_trees_hashes;
        std::size_t proof_size_roots_hash_number =
                1 + 1 + (localization_parameter_array.size() * inter_repetition_parameter);
        std::size_t proof_size_hash_number = proof_size_roots_hash_number + proof_size_path_hash_number;
        double proof_size_hash = double((proof_size_hash_number * hash_ouput_size) / 1024.0 / 8.0);
        total_proof_size_hash += proof_size_hash;

        //std::cout << "proof_size_hash_number is " << proof_size_hash_number << std::endl;
        //std::cout << "secret_vector_tree_hashes is " << secret_vector_tree_hashes*hash_ouput_size/1024/8 << std::endl;
        //std::cout << "h_tree_hashes is " << h_tree_hashes*hash_ouput_size/1024/8 << std::endl;
        //std::cout << "FRI_tree_hashes is " << FRI_trees_hashes*hash_ouput_size/1024/8 << std::endl;
        //std::cout << "proof_size_hash is " << proof_size_hash << std::endl;

        std::size_t proof_size_field_number = 0;
        std::size_t poly_degree_bound_last_round = FRI_degree_bound;
        //std::cout << "FRI_degree_bound 1 " << FRI_degree_bound << std::endl;
        // Note that it is necessary to query (query number) times in all rounds instead of each round
        for (std::size_t i = 0; i < localization_parameter_array.size(); i++) {
            if (i == 0) {
                // v field number
                proof_size_field_number +=
                        (instance + 1) * query_repetition_parameter * (1 << localization_parameter_array[i]);
                //std::cout << "number 1 " << proof_size_field_number << std::endl;
                // h field number
                proof_size_field_number += query_repetition_parameter * (1 << localization_parameter_array[i]);
                //std::cout << "number 2 " << proof_size_field_number << std::endl;
                poly_degree_bound_last_round /= (1 << localization_parameter_array[i]);
                //std::cout << "FRI_degree_bound 2 " << poly_degree_bound_last_round << std::endl;
            } else {
                proof_size_field_number += query_repetition_parameter * (1 << localization_parameter_array[i]) *
                                           inter_repetition_parameter;
                //std::cout << "number 3 " << proof_size_field_number << std::endl;
                poly_degree_bound_last_round /= (1 << localization_parameter_array[i]);
                //std::cout << "FRI_degree_bound 3 " << poly_degree_bound_last_round << std::endl;
            }
        }

        // The final_poly proof size
        std::size_t FRI_field_number_last_round = inter_repetition_parameter * (poly_degree_bound_last_round + 1);
        //std::cout << "number 4 " << FRI_field_number_last_round << std::endl;
        proof_size_field_number += FRI_field_number_last_round;
        double proof_size_field = double((proof_size_field_number * field_size_bits) / 1024.0 / 8.0);
        total_proof_size_field += proof_size_field;
        std::cout << "proof_size_field_number is " << proof_size_field_number << std::endl;
        std::cout << "proof_size_field is " << proof_size_field << std::endl;

        double total_proof_size = proof_size_field + proof_size_hash;
//        libff::print_indent();
//        printf("* total proof size = %f\n", total_proof_size);
    }

    std::cout<<"protocol run correctly!"<<std::endl;
    std::cout<<"prover time, pubpoly time, verifier time, proof size field, proof size hash is "<<std::endl;
    std::cout<< (prover_time)/repeat_num<< '\t'   << (verifier_time)/repeat_num << '\t'
              << (total_proof_size_hash+ total_proof_size_field)/repeat_num << std::endl;
//    std::cout<<"verifier time is "<<std::endl;
//    std::cout<<(verifier_time+pubpoly_time)/repeat_num<<" s"<<std::endl;
//    std::cout<<"pubpoly time is "<<std::endl;
//    std::cout<<pubpoly_time/repeat_num<<" s"<<std::endl;
//    std::cout << "proof size hash is " <<std::endl;
//    std::cout << total_proof_size_hash/repeat_num << std::endl;
//    std::cout << "proof size field is " <<std::endl;
//    std::cout << total_proof_size_field/repeat_num << std::endl;
//    std::cout << "total proof size is " << std::endl;
//    std::cout << (total_proof_size_hash+ total_proof_size_field)/repeat_num << std::endl;
}