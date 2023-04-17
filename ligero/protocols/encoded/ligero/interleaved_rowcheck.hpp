/**@file
 *****************************************************************************
 Interfaces for Interleaved Rowcheck, with oracle target (the target is given
 as an oracle, not a public vector). Tests that messages encoded by
 interleaved RS codes satisfy a given quadratic equation. Based on
 Test-Quadratic-Constraints from Ligero protocol [ACIV17], section 4.3.

 Specifically, checks that oracles of presumed Reed-Solomon codewords,
 encoding messages x, y, and z, satisfy x ⦿ y + a ⦿ z = b, for given a and b
 (where ⦿ denotes pointwise product).
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_PROTOCOLS_ENCODED_LIGERO_INTERLEAVED_ROWCHECK_HPP_
#define ligero_PROTOCOLS_ENCODED_LIGERO_INTERLEAVED_ROWCHECK_HPP_

#include <cstddef>
#include <memory>
#include <vector>

#include "ligero/algebra/polynomials/polynomial.hpp"
#include "ligero/algebra/field_subset/subspace.hpp"
#include "ligero/iop/iop.hpp"

namespace ligero {

template<typename FieldT>
class interleaved_rowcheck_protocol {
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

    std::size_t num_oracles_;
    std::size_t num_oracles_input_;
    std::size_t num_queries_;
    std::size_t zk_queries;

    std::size_t num_interactions_; // sigma in Ligero paper

    const bool make_zk_;
    const field_subset_type field_subset_type_;

    std::vector<std::vector<std::vector<FieldT>>> x_query;
    std::vector<std::vector<std::vector<FieldT>>> y_query;
    std::vector<std::vector<std::vector<FieldT>>> z_query;
    std::vector<std::vector<FieldT>> blinding_query;

public:
    /* Initialization and registration */
    interleaved_rowcheck_protocol(iop_protocol<FieldT> &IOP,
                                  const domain_handle &codeword_domain_handle,
                                  const domain_handle &systematic_domain_handle,
                                  const domain_handle &extended_systematic_domain_handle,
                                  const domain_handle &RS_systematic_domain_handle,
                                  const std::size_t num_oraces_input_,
                                  const std::size_t num_oracles,
                                  const std::size_t num_queries,
                                  const std::size_t num_interactions,
                                  const bool make_zk,
                                  const field_subset_type domain_type);

    /* Proving */
    void calculate_and_submit_responses(const std::vector<std::vector<FieldT>> &Ux_matrix,
                                        const std::vector<std::vector<FieldT>> &Uy_matrix,
                                        const std::vector<std::vector<FieldT>> &Uz_matrix,
                                        const std::vector<std::size_t>& query_set,
                                        const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                        const std::vector<std::vector<FieldT>> &blinding_vectors);

    void calculate_and_submit_responses(const std::vector<std::vector<FieldT>> &Ux_matrix,
                                        const std::vector<std::vector<FieldT>> &Uy_matrix,
                                        const std::vector<std::vector<FieldT>> &Uz_matrix,
                                        const std::vector<std::size_t>& query_set,
                                        const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                        const std::vector<polynomial<FieldT>> public_polys,
                                        const std::vector<std::vector<FieldT>> &blinding_vectors);

    /* Verification */
    bool verifier_predicate(const std::vector<std::vector<FieldT>> &random_vectors,
                            const std::vector<std::vector<FieldT>> &response_mes,
                            std::vector<std::size_t> query_set);
    std::vector<std::vector<FieldT>> response_coefficients;
    std::vector<polynomial<FieldT>> secret_polys;
    std::vector<polynomial<FieldT>> public_polys;
    FieldT target_sum;
};

} // namespace ligero

#include "ligero/protocols/encoded/ligero/interleaved_rowcheck.tcc"

#endif
