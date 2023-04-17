/**@file
 *****************************************************************************
 Utility methods for IOP query position handles
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_IOP_UTILITIES_QUERY_POSITIONS_HPP_
#define ligero_IOP_UTILITIES_QUERY_POSITIONS_HPP_

#include <vector>

#include "ligero/algebra/field_subset/field_subset.hpp"
#include "ligero/iop/iop.hpp"

namespace ligero {

/** TODO: Come up with better name
 *
 */
template<typename FieldT>
std::vector<query_position_handle> query_position_to_queries_for_entire_coset(
    iop_protocol<FieldT> &IOP,
    const query_position_handle &initial_query,
    const field_subset<FieldT> &domain,
    const size_t coset_size);

} // namespace ligero

#include "ligero/iop/utilities/query_positions.tcc"

#endif // ligero_IOP_UTILITIES_QUERY_POSITIONS_HPP_
