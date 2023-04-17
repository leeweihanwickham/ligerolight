/**@file
 *****************************************************************************
 Implementation of exponentiation algorithms.
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_ALGEBRA_EXPONENTIATION_HPP_
#define ligero_ALGEBRA_EXPONENTIATION_HPP_

#include "ligero/algebra/field_subset/field_subset.hpp"
#include "ligero/algebra/field_subset/subgroup.hpp"
#include "ligero/algebra/field_subset/subspace.hpp"

namespace ligero {

template<typename FieldT>
std::vector<FieldT> subset_element_powers(const field_subset<FieldT> &S,
                                          const std::size_t exponent);

template<typename FieldT>
std::vector<FieldT> subspace_element_powers(const affine_subspace<FieldT> &S,
                                            const std::size_t exponent);

template<typename FieldT>
std::vector<FieldT> coset_element_powers(const multiplicative_coset<FieldT> &S,
                                         const std::size_t exponent);

} // namespace ligero

#include "ligero/algebra/exponentiation.tcc"

#endif // ligero_ALGEBRA_EXPONENTIATION_HPP_
