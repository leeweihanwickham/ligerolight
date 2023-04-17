/**@file
 *****************************************************************************
 Ligero protocol for R1CS (from [ACIV17])
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_PROTOCOLS_ENCODED_LIGERO_LIGERO_HPP_
#define ligero_PROTOCOLS_ENCODED_LIGERO_LIGERO_HPP_

#include <cstddef>
#include <memory>
#include <vector>
#include <algorithm>

#include "ligero/iop/iop.hpp"
#include "ligero/relations/r1cs.hpp"

#include <libff/algebra/field_utils/field_utils.hpp>
#include "ligero/algebra/field_subset/subgroup.hpp"
#include "ligero/algebra/polynomials/polynomial.hpp"
#include "ligero/algebra/field_subset/subspace.hpp"
#include "ligero/protocols/encoded/ligero/interleaved_lincheck_ot.hpp"
#include "ligero/protocols/encoded/ligero/interleaved_rowcheck.hpp"
#include "ligero/protocols/encoded/ligero/interactive_RS_test.hpp"
#include "ligero/protocols/ldt/fri/fri_ldt.hpp"
/* commitment */
#include "ligero/bcs/hashing/blake2b.hpp"
#include "ligero/bcs/hashing/dummy_algebraic_hash.hpp"
#include "ligero/bcs/hashing/hashing.hpp"
#include "ligero/bcs/merkle_tree.hpp"
#include "ligero/bcs/Newmerkle.hpp"
#include "ligero/protocols/fft_circuit_GKR.hpp"
namespace ligero {

struct encoded_ligero_parameters {
    std::size_t num_interaction_phase_repetitions_;
    std::size_t num_query_phase_repetitions_;

    std::size_t IPA_num_interaction_repetitions_;
    std::size_t IPA_num_query_repetitions_;

    bool make_zk_;
    field_subset_type field_subset_type_;

    std::size_t matrix_width_;
    std::size_t matrix_height_;
    std::size_t num_oracles_input_;
    std::size_t num_oracles_vectors_;
    std::vector<std::size_t> localization_parameter_array;
    std::size_t col_systematic_domain_dim_;
};

template<typename FieldT>
class interleaved_r1cs_protocol {
protected:

    std::vector<std::vector<FieldT>> matrix_w_commit;
    std::vector<std::vector<FieldT>> matrix_x_commit;
    std::vector<std::vector<FieldT>> matrix_y_commit;
    std::vector<std::vector<FieldT>> matrix_z_commit;
    std::vector<std::vector<FieldT>> matrix_add_commit;
    std::vector<std::vector<FieldT>> matrix_Uxyz_Rowcheck_blinding_commit;

    iop_protocol<FieldT> &IOP_;
    r1cs_constraint_system<FieldT> constraint_system_;
    encoded_ligero_parameters parameters_;

    domain_handle codeword_domain_handle_;
    domain_handle systematic_domain_handle_;
    domain_handle extended_systematic_domain_handle_;
    domain_handle RS_col_systematic_domain_handle_;
    domain_handle RS_col_codeword_domain_handle_;

    field_subset<FieldT> codeword_domain_;
    field_subset<FieldT> systematic_domain_;
    field_subset<FieldT> extended_systematic_domain_;
    field_subset<FieldT> RS_col_systematic_domain_;
    field_subset<FieldT> RS_col_codeword_domain_;

    std::size_t codeword_domain_size_;
    std::size_t systematic_domain_size_;
    std::size_t extended_systematic_domain_size_;
    std::size_t RS_col_systematic_domain_size_;
    std::size_t RS_col_codeword_domain_size_;

    std::shared_ptr<interleaved_lincheck_ot_protocol<FieldT> > lincheck_A_;
    std::shared_ptr<interleaved_lincheck_ot_protocol<FieldT> > lincheck_B_;
    std::shared_ptr<interleaved_lincheck_ot_protocol<FieldT> > lincheck_C_;
    std::shared_ptr<interleaved_rowcheck_protocol<FieldT> > rowcheck_;
    std::shared_ptr<interleaved_RScode<FieldT>> IRScheck;
    std::shared_ptr<interleaved_RScode<FieldT>> IRScheck_w_;
    std::shared_ptr<interleaved_RScode<FieldT>> IRScheck_x_;
    std::shared_ptr<interleaved_RScode<FieldT>> IRScheck_y_;
    std::shared_ptr<interleaved_RScode<FieldT>> IRScheck_z_;
    std::shared_ptr<Inner_product_verifier<FieldT>> IPA_verifier_;
    std::shared_ptr<Inner_product_prover<FieldT>> IPA_prover_;
    std::shared_ptr<merkle<FieldT>> merkleTree;

    //mark
    std::vector<fft_GKR<FieldT>> fft_GKRs;

    std::size_t matrix_height_;
    std::size_t matrix_width_;

    std::size_t num_oracles_input_;
    std::size_t num_oracles_vectors_;
    std::size_t num_queries_;
    std::size_t num_interactions_;
    //std::size_t Allmerklelenth=0;
    bool make_zk_;
    field_subset_type field_subset_type_;
    const std::size_t parallel_opt_;
    std::size_t encoding_independence_;

    std::vector<std::vector<FieldT>> lincheck_blinding_vector;
    std::vector<std::vector<FieldT>> rowcheck_blinding_vector;

    naive_sparse_matrix<FieldT> A_matrix_;
    naive_sparse_matrix<FieldT> B_matrix_;
    naive_sparse_matrix<FieldT> C_matrix_;
//  未按列编码 只需要一轮即可
    std::vector<std::vector<FieldT>> Uw_matrix;
    std::vector<std::vector<FieldT>> Ux_matrix;
    std::vector<std::vector<FieldT>> Uy_matrix;
    std::vector<std::vector<FieldT>> Uz_matrix;
    std::vector<std::vector<FieldT>> additional_matrix;
//blinding矩阵需要多轮
    std::vector<std::vector<FieldT>> rowcheck_blinding_matrix;
    /* commitment paramaters */
    std::vector<std::vector<FieldT>> value_matrix_for_commitment;

    // 按列编码后的矩阵
//    std::vector<std::vector<FieldT>> matrix_w_commit;
//    std::vector<std::vector<FieldT>> matrix_x_commit;
//    std::vector<std::vector<FieldT>> matrix_y_commit;
//    std::vector<std::vector<FieldT>> matrix_z_commit;
//    std::vector<std::vector<FieldT>> matrix_add_commit;
//    std::vector<std::vector<FieldT>> matrix_Uxyz_Rowcheck_blinding_commit;

    std::vector<std::vector<FieldT>> blinding_matrix_for_commitment;

    std::vector<merkle_tree_set_membership_proof<std::string>> commitments;

    std::vector<std::size_t> query_set;
    std::vector<std::size_t> IPA_query_set;
    std::vector<std::size_t> query_set_sort_by_sequence;
    std::vector<std::size_t> query_col_set;
    std::vector<std::vector<std::vector<FieldT>>> query_value;
//  重构的merkle
    std::vector<merkleTreeParameter> lincheck_blinding_commit;
    std::vector<std::vector<FieldT>> IRS_random_combinations;// 每轮r：IRS_random_combinations[i]
    std::vector<std::vector<FieldT>> lincheck_ABC_randomcombinations;
    std::vector<std::vector<FieldT>> rowcheck_randomcombinations;
    const std::size_t IPA_RS_extra_dimensions_;
public:
    /* Initialization and registration */
    interleaved_r1cs_protocol(iop_protocol<FieldT> &IOP,
                              const domain_handle &codeword_domain_handle,
                              const domain_handle &systematic_domain_handle,
                              const domain_handle &extended_systematic_domain_handle,
                              const domain_handle &RS_col_systematic_domain_handle,
                              const domain_handle &RS_col_codeword_domain_handle,
                              const r1cs_constraint_system<FieldT> &constraint_system,
                              const encoded_ligero_parameters &parameters,
                              const std::size_t parallel_opt,
                              const std::size_t IPA_RS_extra_dimensions);


    void register_queries();
    //void register_col_queries();
    /* Proving */
    float submit_witness_oracles(const std::vector<FieldT> &primary_input,
                                const std::vector<FieldT> &auxiliary_input,
                                const std::map<std::size_t,FieldT> &public_input);
    void submit_blinding_vector_oracles();

    float inner_product_prover(const encoded_ligero_parameters &parameters);

    bool inner_product_verifier(const encoded_ligero_parameters &parameters);

//    bool consistency_check();

    float calculate_and_submit_proof(const std::map<std::size_t,FieldT> &public_input);

    /* Verification */
    bool verifier_predicate(const std::map<std::size_t,FieldT> &public_input);

protected:
    std::vector<FieldT> submit_zero_sum_blinding_vector();
    std::vector<FieldT> submit_zero_blinding_vector();
};

} // namespace ligero

#include "ligero/protocols/encoded/ligero/ligero.tcc"

#endif // ligero_PROTOCOLS_ENCODED_LIGERO_LIGERO_HPP_

