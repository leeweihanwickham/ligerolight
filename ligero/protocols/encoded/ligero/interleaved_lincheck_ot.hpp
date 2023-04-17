/**@file
 *****************************************************************************
 Interfaces for Interleaved Lincheck, with oracle target (the target is given
 as an oracle, not a public vector). Tests that messages encoded by interleaved
 RS codes satisfy a given linear relation.
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_PROTOCOLS_ENCODED_LIGERO_INTERLEAVED_LINCHECK_OT_HPP_
#define ligero_PROTOCOLS_ENCODED_LIGERO_INTERLEAVED_LINCHECK_OT_HPP_

#include <cstddef>
#include <memory>
#include <vector>

#include <libff/algebra/field_utils/field_utils.hpp>
#include "ligero/algebra/field_subset/subgroup.hpp"
#include "ligero/algebra/polynomials/polynomial.hpp"
#include "ligero/algebra/field_subset/subspace.hpp"
#include "ligero/iop/iop.hpp"
#include "ligero/relations/r1cs.hpp"

namespace ligero {

template<typename FieldT>
class interleaved_lincheck_ot_protocol {
protected:
    iop_protocol<FieldT> &IOP_;
    domain_handle codeword_domain_handle_;
    domain_handle systematic_domain_handle_;
    domain_handle RS_col_systematic_domain_handle_;

    field_subset<FieldT> codeword_domain_;
    field_subset<FieldT> systematic_domain_;
    field_subset<FieldT> extended_systematic_domain_;
    field_subset<FieldT> RS_col_systematic_domain_;

    std::size_t codeword_domain_size_;
    std::size_t systematic_domain_size_;
    std::size_t RS_col_systematic_domain_size_;
    std::size_t response_size_;

    const std::size_t num_oracles_input_;
    const std::size_t num_oracles_target_;
    const std::size_t num_queries_;
    std::size_t zk_queries;
    const std::size_t parallel_opt_;
    const std::size_t num_interactions_; // sigma in Ligero paper

    const bool make_zk_;
    const field_subset_type field_subset_type_;

    const naive_sparse_matrix<FieldT> constraint_matrix_;

    std::vector<verifier_random_message_handle> random_linear_combination_handles_;
    std::vector<prover_message_handle> response_handles_;

    std::vector<oracle_handle_ptr> input_vector_row_oracle_handles_;
    std::vector<oracle_handle_ptr> target_vector_row_oracle_handles_;
    std::vector<oracle_handle_ptr> blinding_vector_row_oracle_handles_;

    std::vector<random_query_position_handle> query_position_handles_;
    std::vector<std::vector<query_handle>> input_vector_row_oracles_by_position_handles_;
    std::vector<std::vector<query_handle>> target_vector_row_oracles_by_position_handles_;
    std::vector<std::vector<query_handle>> blinding_vector_row_oracles_by_position_handles_;

public:
    /* Initialization and registration */
    interleaved_lincheck_ot_protocol(iop_protocol<FieldT> &IOP,
                                     const domain_handle &codeword_domain_handle,
                                     const domain_handle &systematic_domain_handle,
                                     const domain_handle &extended_systematic_domain_handle,
                                     const domain_handle &RS_col_systematic_domain_handle,
                                     const std::size_t num_oracles_input,
                                     const std::size_t num_oracles_target,
                                     const std::size_t num_queries,
                                     const std::size_t num_interactions,
                                     const bool make_zk,
                                     const field_subset_type domain_type,
                                     const naive_sparse_matrix<FieldT> constraint_matrix,
                                     const std::size_t parallel_opt);

    void calculate_and_submit_responses(const std::vector<std::vector<FieldT>>&addition_matrix,
                                        const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                        const std::vector<std::vector<FieldT>> &blinding_vectors,
                                        const std::vector<std::size_t>& query_set,
                                        const std::vector<std::vector<FieldT>> &input_matrix,
                                        const std::vector<std::vector<FieldT>>&target_matrix);


    /* Verification */
//    bool verifier_predicate();
    bool verifier_predicate(const std::vector<FieldT> &supplementary_input,
                            std::vector<std::vector<FieldT>> random_linear_combinations,
                            const std::vector<std::vector<FieldT>> &response_mes,
                            const std::vector<std::size_t>& query_set);
    std::vector<std::vector<FieldT>> supplementary_input_vectors;
    std::vector<std::vector<FieldT>> response_polys_coefficients;
    std::vector<std::vector<std::vector<FieldT>>> Ur_interactions;
    std::vector<std::vector<std::vector<FieldT>>> input_query_col;
    std::vector<std::vector<std::vector<FieldT>>> target_query_col;
    std::vector<polynomial<FieldT>> secret_polys;
    std::vector<polynomial<FieldT>> public_polys;
    FieldT target_sum;
};

} // namespace ligero

#include "ligero/protocols/encoded/ligero/interleaved_lincheck_ot.tcc"

#endif // ligero_PROTOCOLS_ENCODED_LIGERO_INTERLEAVED_LINCHECK_OT_HPP_
