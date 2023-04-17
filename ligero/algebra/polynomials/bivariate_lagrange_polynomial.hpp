/**@file
 *****************************************************************************
 Class for bivariate lagrange polynomials.
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_ALGEBRA_POLYNOMIALS_BIVARIATE_LAGRANGE_POLYNOMIAL_HPP_
#define ligero_ALGEBRA_POLYNOMIALS_BIVARIATE_LAGRANGE_POLYNOMIAL_HPP_

#include <cstddef>
#include <vector>

#include "ligero/algebra/field_subset/field_subset.hpp"
#include "ligero/algebra/polynomials/lagrange_polynomial.hpp"
#include "ligero/algebra/polynomials/vanishing_polynomial.hpp"

namespace ligero {

/** Bivariate lagrange polynomials in this codebase refer to succinct polynomials
 *  that are low degree extensions of the unnormalized lagrange basis, evaluated at a given position.
 *
 *  The unnormalized lagrange polynomial basis for a domain S is the set:
 *      { Z_S(X) / (X - s) : s \in S }
 *  This basis has the property that a basis element is 0 at all but one element in S.
 *  The low degree extension of the above is the bivariate polynomial lagrange polynomial:
 *      f(X, Y) = (Z_S(X) - Z_S(Y)) / (X - Y)
 *  This low degree extension identifies elements in S with Y.
 *  We can see that f is degree |S| - 1 in Y,
 *  as (X - Y) divides the numerator from definition of vanishing polynomials.
 */
template<typename FieldT>
class bivariate_lagrange_polynomial {
protected:
    field_subset<FieldT> S_;
    vanishing_polynomial<FieldT> Z_S_;
public:
    explicit bivariate_lagrange_polynomial() {};
    explicit bivariate_lagrange_polynomial(const field_subset<FieldT> &S);

    FieldT evaluation_at_point(const FieldT &x, const FieldT &y) const;
    lagrange_polynomial<FieldT> fix_x(const FieldT &c) const;
    std::vector<FieldT> evaluations_over_field_subset(const FieldT &x, const field_subset<FieldT> &S) const;
};

} // namespace ligero

#include "ligero/algebra/polynomials/bivariate_lagrange_polynomial.tcc"

#endif // ligero_ALGEBRA_POLYNOMIALS_BIVARIATE_LAGRANGE_POLYNOMIAL_HPP_
