/**@file
 *****************************************************************************
 Computation of Lagrange coefficients for additive subspaces.
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_ALGEBRA_LAGRANGE_HPP_
#define ligero_ALGEBRA_LAGRANGE_HPP_

#include <vector>

#include "ligero/algebra/field_subset/field_subset.hpp"
#include "ligero/algebra/field_subset/subspace.hpp"
#include "ligero/algebra/polynomials/polynomial.hpp"
#include "ligero/algebra/polynomials/vanishing_polynomial.hpp"

namespace ligero {

template<typename FieldT>
class lagrange_cache {
private:
    void construct_internal(const affine_subspace<FieldT> &domain);
    void construct_internal(const multiplicative_coset<FieldT> &domain);
    std::vector<FieldT> subspace_coefficients_for(const FieldT &interpolation_point);
    std::vector<FieldT> coset_coefficients_for(const FieldT &interpolation_point);
protected:
    const field_subset<FieldT> domain_;
    const vanishing_polynomial<FieldT> vp_;
    // used in additive case
    FieldT c_;
    // used in multiplicative case
    std::vector<FieldT> v_inv_;

    // Cache last evaluation. Does this need to cache more than that?
    const bool cache_evaluations_;
    bool called_;
    FieldT last_interpolation_point_;
    std::vector<FieldT> last_interpolation_result_;

    // Optimization if the interpolation domain doesn't intersect with the base domain
    const bool interpolation_domain_intersects_domain_;
public:
    lagrange_cache(const field_subset<FieldT> &domain,
                   const bool cache_evaluations = false,
                   const bool interpolation_domain_intersects_domain = false);
    lagrange_cache(const affine_subspace<FieldT> &domain,
                   const bool cache_evaluations = false,
                   const bool interpolation_domain_intersects_domain = false);
    lagrange_cache(const multiplicative_coset<FieldT> &domain,
                   const bool cache_evaluations = false,
                   const bool interpolation_domain_intersects_domain = false);
    std::vector<FieldT> coefficients_for(const FieldT &interpolation_point);
};

template<typename FieldT>
std::vector<FieldT> lagrange_coefficients(const field_subset<FieldT> &domain,
                                          const FieldT &interpolation_point);

template<typename FieldT>
std::vector<FieldT> lagrange_coefficients(const affine_subspace<FieldT> &domain,
                                          const FieldT &interpolation_point);

template<typename FieldT>
std::vector<FieldT> lagrange_coefficients(const multiplicative_coset<FieldT> &domain,
                                          const FieldT &interpolation_point);

} // namespace ligero

#include "ligero/algebra/lagrange.tcc"

#endif // ligero_ALGEBRA_LAGRANGE_HPP_
