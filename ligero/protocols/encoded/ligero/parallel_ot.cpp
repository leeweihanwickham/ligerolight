#include <vector>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/fields/prime_base/fields_64.hpp>
#include "ligero/algebra/field_subset/field_subset.hpp"
#include "ligero/algebra/fft.hpp"

//field 128

//void parallel_ot(std::vector<std::vector<libff::alt_bn128_Fr>> &randomized_matrix_row_vector_evals,
//              std::vector<std::vector<libff::alt_bn128_Fr>> &randomized_matrix_row_vectors,
//              ligero::field_subset<libff::alt_bn128_Fr> &systematic_domain_,
//              ligero::field_subset<libff::alt_bn128_Fr> &codeword_domain_,
//              int index,
//              int limit){
//    for(int i=index;i<limit+index;i++){
//        std::vector<libff::alt_bn128_Fr> row_vector_poly_coeffs= IFFT_over_field_subset(randomized_matrix_row_vectors[i],systematic_domain_);
//        std::vector<libff::alt_bn128_Fr> row_vector_poly_evals= FFT_over_field_subset(row_vector_poly_coeffs,codeword_domain_);
//        randomized_matrix_row_vector_evals[i]=std::move(row_vector_poly_evals);
//    }
//}

// field 64
void parallel_ot(std::vector<std::vector<libff::Fields_64>> &randomized_matrix_row_vector_evals,
                 std::vector<std::vector<libff::Fields_64>> &randomized_matrix_row_vectors,
                 ligero::field_subset<libff::Fields_64> &systematic_domain_,
                 ligero::field_subset<libff::Fields_64> &codeword_domain_,
                 int index,
                 int limit){
    for(int i=index;i<limit+index;i++){
        std::vector<libff::Fields_64> row_vector_poly_coeffs= IFFT_over_field_subset(randomized_matrix_row_vectors[i],systematic_domain_);
        std::vector<libff::Fields_64> row_vector_poly_evals= FFT_over_field_subset(row_vector_poly_coeffs,codeword_domain_);
        randomized_matrix_row_vector_evals[i]=std::move(row_vector_poly_evals);
    }
}





