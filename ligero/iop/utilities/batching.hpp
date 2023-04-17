/**@file
 *****************************************************************************
 Utilities for applying the same IOP operation many times
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_IOP_UTILITIES_BATCHING_HPP_
#define ligero_IOP_UTILITIES_BATCHING_HPP_

#include <vector>
#include "ligero/iop/iop.hpp"

namespace ligero {

template<typename FieldT>
std::vector<std::shared_ptr<std::vector<FieldT>>> get_all_oracle_evaluations(
    iop_protocol<FieldT> &IOP,
    const std::vector<oracle_handle_ptr> poly_handles);

template<typename FieldT>
std::vector<prover_message_handle> register_n_prover_messages(
    iop_protocol<FieldT> &IOP,
    const size_t num_messages,
    const size_t message_size);

template<typename FieldT>
std::vector<verifier_random_message_handle> register_n_verifier_messages(
    iop_protocol<FieldT> &IOP,
    const size_t num_messages,
    const size_t message_size);

template<typename FieldT>
std::vector<oracle_handle_ptr> register_n_oracles(
    iop_protocol<FieldT> &IOP,
    const size_t num_oracles,
    const domain_handle domain,
    const size_t degree_bound,
    const bool make_zk);

std::vector<oracle_handle_ptr> virtual_oracle_handles_to_handle_ptrs(
    const std::vector<virtual_oracle_handle> handles);

template<typename FieldT>
std::vector<query_handle> register_queries_for_same_pos(
    iop_protocol<FieldT> &IOP,
    const std::vector<oracle_handle_ptr> oracle_handles,
    const query_position_handle query_position);

} // namespace ligero

#include "ligero/iop/utilities/batching.tcc"

#endif // ligero_IOP_UTILITIES_BATCHING_HPP_
