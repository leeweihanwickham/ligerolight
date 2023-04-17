/**@file
 *****************************************************************************
 Utility functions for subspace bases.
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_ALGEBRA_BASIS_UTILS_HPP_
#define ligero_ALGEBRA_BASIS_UTILS_HPP_

#include <cstddef>
#include <vector>
#include <libff/algebra/field_utils/field_utils.hpp>
#include "ligero/algebra/polynomials/polynomial.hpp"

namespace ligero {

/** Creates a basis of monomials, of the specified dimension.
 *  This method returns the basis x^i, x^{i + 1}, ..., x^{i + d - 1},
 *  where i is the smallest exponent, and d is the dimension.
 */
template<typename FieldT>
std::vector<FieldT> monomial_basis(const size_t dimension, const size_t smallest_exponent);

/** Returns the evaluations of the provided polynomial over the basis. */
template<typename FieldT>
std::vector<FieldT> transform_basis_by_polynomial(
    const std::shared_ptr<polynomial_base<FieldT>> transform, const std::vector<FieldT> basis);

} // namespace ligero

#include "ligero/algebra/field_subset/basis_utils.tcc"

#endif // ligero_ALGEBRA_BASIS_UTILS_HPP_
