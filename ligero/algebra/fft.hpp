/**@file
 *****************************************************************************
 Implementation of Gao-Mateer for the additive FFT/IFFT,
 implementation of Nlog(d) Cooley-Tukey for the multiplicative FFT,
 and wrappers to libfqfft for the multiplicative IFFT.
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_ALGEBRA_FFT_HPP_
#define ligero_ALGEBRA_FFT_HPP_

#include <vector>

#include "ligero/algebra/field_subset/field_subset.hpp"
#include "ligero/algebra/field_subset/subspace.hpp"
#include "ligero/algebra/field_subset/subgroup.hpp"
#include <libff/common/utils.hpp>

namespace ligero {

/* Performs naive computation of the polynomial evaluation
   problem. Mostly useful for testing. */
template<typename FieldT>
std::vector<FieldT> naive_FFT(const std::vector<FieldT> &poly_coeffs,
                              const field_subset<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                 const affine_subspace<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals,
                                  const affine_subspace<FieldT> &domain);

/* Calls additive_FFT but adds trace data */
template<typename FieldT>
std::vector<FieldT> additive_FFT_wrapper(const std::vector<FieldT> &v,
                                         const affine_subspace<FieldT> &H);

/* Calls additive_IFFT but adds trace data */
template<typename FieldT>
std::vector<FieldT> additive_IFFT_wrapper(const std::vector<FieldT> &v,
                                          const affine_subspace<FieldT> &H);

template<typename FieldT>
std::vector<FieldT> multiplicative_FFT(const std::vector<FieldT> &poly_coeffs,
                                       const multiplicative_coset<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> multiplicative_IFFT(const std::vector<FieldT> &evals,
                                        const multiplicative_coset<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> multiplicative_FFT_wrapper(const std::vector<FieldT> &v,
                                               const multiplicative_coset<FieldT> &H);

template<typename FieldT>
std::vector<FieldT> multiplicative_IFFT_wrapper(const std::vector<FieldT> &v,
                                                const multiplicative_coset<FieldT> &H);

template<typename FieldT>
std::vector<FieldT> FFT_over_field_subset(const std::vector<typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type> coeffs,
                                            field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> FFT_over_field_subset(const std::vector<typename libff::enable_if<libff::is_additive<FieldT>::value, FieldT>::type> coeffs,
                                            field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> IFFT_over_field_subset(const std::vector<typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type> evals,
                                             field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> IFFT_over_field_subset(const std::vector<typename libff::enable_if<libff::is_additive<FieldT>::value, FieldT>::type> evals,
                                             field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> IFFT_of_known_degree_over_field_subset(
    const std::vector<typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type> evals,
    size_t degree_bound,
    field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> IFFT_of_known_degree_over_field_subset(
    const std::vector<typename libff::enable_if<libff::is_additive<FieldT>::value, FieldT>::type> evals,
    size_t degree_bound,
    field_subset<FieldT> domain);

} // namespace ligero

#include "ligero/algebra/fft.tcc"

#endif // ligero_ALGEBRA_FFT_HPP_
