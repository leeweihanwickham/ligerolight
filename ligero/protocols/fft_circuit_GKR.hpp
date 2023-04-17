#ifndef fft_circuit_GKR
#define fft_circuit_GKR

#include <cstdio>
#include <vector>
#include <chrono>

#include <cstddef>

#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>
#include <libfqfft/evaluation_domain/domains/basic_radix2_domain_aux.hpp>

#include <libff/common/profiling.hpp>
#include <libff/algebra/field_utils/field_utils.hpp>
#include "ligero/algebra/utils.hpp"

#include "ligero/algebra/field_subset/field_subset.hpp"
#include "ligero/algebra/field_subset/subspace.hpp"
#include "ligero/algebra/field_subset/subgroup.hpp"
#include <libff/common/utils.hpp>


#include "ligero/protocols/fft_circuit_GKR_polynomial.hpp"

namespace ligero{

template<typename FieldT>
class circuit
{
public:
    std::vector<FieldT*> circuit_val;
    std::vector<int> size;
    ~circuit(){
        for(size_t i = 0; i < circuit_val.size(); i++)
        {
            delete circuit_val[i];
        }
    }
};


template<typename FieldT>
class fft_GKR{
protected:
    void init_array(int max_bit_length, bool prover = true);

    bool destruction_array(bool prover = true);

public:
    field_subset<FieldT> coset;

    /* fft parameters */
    //circuit<FieldT> C;   // circuit
    //std::vector<FieldT> poly_coeffs;
    bool invert;  // fft:invert = 0, ifft:invert = 1;
    FieldT rou, inv_rou; // the root/inv_root of coset

    FieldT inv_n;
    int lg_size, poly_dimension;
    size_t poly_size;
    FieldT *eval_points;
    //std::vector<FieldT> result; // fft/ifft result

    static const size_t num_poly = 64;

    /* GKR parameters */
    std::vector<FieldT> alpha_list, beta_list;
    FieldT alpha_beta_sum0;
    std::vector<FieldT> r_0, one_minus_r_0, r_1, one_minus_r_1;

    std::vector<FieldT> v_v_list, v_u_list;
    std::vector<std::vector<FieldT>> r_u_list, r_v_list;


    std::vector<quadratic_poly<FieldT>> addition_layer_poly, mult_layer_poly;
    std::vector<std::vector<quadratic_poly<FieldT>>> fft_layer_phase1_poly, fft_layer_phase2_poly;

    //temp arrays
    FieldT *beta_g_r0_fhalf, *beta_g_r0_shalf, *beta_g_r1_fhalf, *beta_g_r1_shalf, *beta_u_fhalf, *beta_u_shalf;
    linear_poly<FieldT> *add_mult_sum, *V_mult_add, *addV_array; // TODO: replace these with algerbra/polynomial


    /* evaluate the protocol */
    size_t proof_size;
    float v_time, p_time;
    float p_fft_time;

    fft_GKR(const field_subset<FieldT>coset,
            //const std::vector<FieldT> poly_coeffs_,
            bool invert_ = 0);

    /*
    ~fft_GKR(){
        if(this->eval_points != NULL)
            delete[] this->eval_points;
    }
    */
    /* initialize  functions */
    void build_circuit(circuit<FieldT> &C, const std::vector<FieldT> &poly_coeffs,std::vector<FieldT> &result);


    /* GKR stages */

    FieldT V_output(const FieldT *output_raw, int r_0_size, int output_size);

    bool addition_layer(FieldT *c_val, int num_poly, FieldT &alpha_beta_sum);

    bool mult_layer(FieldT *c_val, int num_poly, FieldT &alpha_beta_sum);

    bool intermediate_layer(FieldT &alpha_beta_sum);

    quadratic_poly<FieldT> sumcheck_phase1_update(FieldT previous_random, int current_bit, int total_uv);
    quadratic_poly<FieldT> sumcheck_phase2_update(FieldT previous_random, int current_bit, int total_uv);

    bool extension_gkr();

    bool ifft_gkr(circuit<FieldT> &C, FieldT &alpha_beta_sum);

    int engage_gkr(const std::vector<FieldT> &poly_coeffs,std::vector<FieldT> &result);

    bool prover_compute(const std::vector<FieldT> &poly_coeffs,std::vector<FieldT> &result);

    /* verifier computation */
    bool addition_layer_verify(FieldT &alpha_beta_sum);

    bool mult_layer_verify(FieldT &alpha_beta_sum);

    bool ifft_gkr_verify(FieldT &alpha_beta_sum);

    bool verifier_predicate();


    /* auxiliary functions */
    void refresh_randomness(FieldT *r, FieldT *one_minus_r, int size);

    void refresh_randomness(std::vector<FieldT> &r, std::vector<FieldT> &one_minus_r, int size);

};

template<typename FieldT>
void random_field_vector(std::vector<FieldT> &v, size_t n);

}

#include <ligero/protocols/fft_circuit_GKR.tcc>
#endif