/**@file
 *****************************************************************************
 Specialized IOP that implements the BCS16 transformation.
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_SNARK_COMMON_COMMON_BCS_PARARMETERS_HPP_
#define ligero_SNARK_COMMON_COMMON_BCS_PARARMETERS_HPP_

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <map>
#include <vector>

#include "ligero/bcs/bcs_common.hpp"
#include "ligero/bcs/hashing/hash_enum.hpp"
#include "ligero/bcs/hashing/hashing.hpp"
#include "ligero/bcs/hashing/blake2b.hpp"
#include "ligero/bcs/pow.hpp"

namespace ligero {

template<typename FieldT, typename MT_root_hash>
bcs_transformation_parameters<FieldT, MT_root_hash> default_bcs_params(
    const bcs_hash_type hash_type, const size_t security_parameter);

} // namespace ligero

#include "ligero/bcs/common_bcs_parameters.tcc"

#endif // ligero_SNARK_COMMON_COMMON_BCS_PARARMETERS_HPP_
