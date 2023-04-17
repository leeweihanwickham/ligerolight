//
// Created by wx on 22-11-26.
//

#ifndef LIGERO_FFT_CIRCUIT_GKR_TCC
#define LIGERO_FFT_CIRCUIT_GKR_TCC
namespace ligero{

template<typename FieldT>
fft_GKR<FieldT>::fft_GKR(const field_subset<FieldT> coset_,
                         //const std::vector<FieldT> poly_coeffs_ ,
                         bool invert_):
        coset(coset_),invert(invert_)
{
    v_time = 0;
    p_fft_time = 0;
    p_time = 0;
    proof_size = 0;
    eval_points = NULL;
    //build_circuit();

}

template<typename FieldT>
void fft_GKR<FieldT>::build_circuit(circuit<FieldT> &C,const std::vector<FieldT> &poly_coeffs,std::vector<FieldT> &result)
{
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    if(invert)
        assert(poly_coeffs.size() == coset.num_elements());
    else
        assert(poly_coeffs.size() <= coset.num_elements());
    const size_t n = coset.num_elements();
    lg_size = libff::log2(n);
    FieldT shift = coset.shift();
    result.resize(n);

    std::vector<FieldT> poly_coeffs_(poly_coeffs);
    if (shift != FieldT::one() && (invert == 0))
    {
        libfqfft::_multiply_by_coset<FieldT>(poly_coeffs_, shift);
    }
    poly_coeffs_.resize(n, FieldT::zero());


    // 仅用来test
    poly_dimension = libff::log2(poly_coeffs.size());
    poly_size = poly_coeffs.size();
    const size_t duplicity_of_initial_elems = 1ull << (lg_size - poly_dimension);

    poly_coeffs_.resize(n, FieldT::zero());
    if (duplicity_of_initial_elems > 1)
    {
        for(size_t k = 1; k < duplicity_of_initial_elems; k++)
        {
            for(size_t j = 0; j < poly_size; j++)
            {
                poly_coeffs_[k * poly_size + j] = poly_coeffs_[j];
            }
        }
        poly_size = n;
    }

    int i = lg_size-1;
    // FFT/IFFT input
    C.circuit_val.push_back(new FieldT[1 << (i + 1)]);
    for(int j = 0; j < (1 << (i+1)); ++j)
    {
        C.circuit_val[0][j] = poly_coeffs_[j];
    }
    C.size.push_back(1 << (i + 1));

    //ifft
    rou = coset.generator();
    assert(libff::power(rou, n) == FieldT::one());
    inv_rou = rou.inverse();
    inv_n = FieldT(n);
    inv_n.invert();
    assert(FieldT::one() == (inv_n * FieldT(n)));

    std::vector<FieldT> x_arr;
    x_arr.resize(1 << lg_size, FieldT::zero());
    std::vector<FieldT> rot_mul;
    rot_mul.resize(62,0);
    rot_mul[0] = invert? inv_rou : rou;
    for(int i = 1; i < 62; ++i)
    {
        rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1];
    }
    int starting_depth = 0;
    for(int dep = poly_dimension - 1; dep >= 0; --dep)
    {
        int blk_size = 1 << (lg_size - dep);
        int half_blk_size = blk_size >> 1;
        int cur = starting_depth + (poly_dimension - dep);
        int pre = cur - 1;
        C.circuit_val.push_back(new FieldT[1 << lg_size]);
        C.size.push_back(1 << lg_size);
        FieldT x = FieldT(1);
        {
            x_arr[0] = FieldT(1);
            for(int j = 1; j < blk_size; ++j)
                x_arr[j] = x_arr[j - 1] * rot_mul[dep];
            for(int k = 0; k < blk_size / 2; ++k)
            {
                int double_k = (k) & (half_blk_size - 1);
                for(int j = 0; j < (1 << dep); ++j)
                {
                    auto l_value = C.circuit_val[pre][double_k << (dep + 1) | j], r_value = x_arr[k] * C.circuit_val[pre][double_k << (dep + 1) | (1 << dep) | j];
                    C.circuit_val[cur][k << dep | j] = l_value + r_value;
                    C.circuit_val[cur][(k + blk_size / 2) << dep | j] = l_value - r_value;
                }
            }
        }
    }
    int last_layer_id = C.circuit_val.size();
    C.circuit_val.push_back(new FieldT[1 << lg_size]);
    C.size.push_back(1 << lg_size);
    if(invert)
    {
        for(int i = 0; i < (1 << lg_size); ++i)
        {
            C.circuit_val[last_layer_id][i] = C.circuit_val[last_layer_id - 1][i] * inv_n;
            result[i] = C.circuit_val[last_layer_id][i];
        }
    }
    else
    {
        for(int i = 0; i < (1 << lg_size); ++i)
        {
            C.circuit_val[last_layer_id][i] = C.circuit_val[last_layer_id - 1][i];
            result[i] = C.circuit_val[last_layer_id][i];
        }
    }

    if (shift != FieldT::one() && invert)
    {
        libfqfft::_multiply_by_coset<FieldT>(result, shift.inverse());
    }
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    p_fft_time = time_span.count();
    //poly eval assuming 64 repetitions
    //directly evaluate the polynomial for practical purpose
    //since the number of repetition is only 40
    FieldT eval_point;
    last_layer_id = C.circuit_val.size();
    C.circuit_val.push_back(new FieldT[64 * (1 << lg_size)]);
    C.size.push_back(64 << lg_size);
    eval_points = new FieldT[64];
    for(int i = 0; i < 64; ++i)
    {
        eval_point = FieldT::random_element();
        eval_points[i] = eval_point;
        FieldT x = FieldT(1);
        for(int j = 0; j < (1 << lg_size); ++j)
        {
            C.circuit_val[last_layer_id][j + (i << lg_size)] = C.circuit_val[last_layer_id - 1][j] * x;
            x = x * eval_point;
        }
    }
    last_layer_id = C.circuit_val.size();
    C.circuit_val.push_back(new FieldT[64]);
    C.size.push_back(64);
    for(int i = 0; i < 64; ++i)
    {
        C.circuit_val[last_layer_id][i] = FieldT(0);
        for(int j = 0; j < (1 << lg_size); ++j)
        {
            C.circuit_val[last_layer_id][i] = C.circuit_val[last_layer_id][i] + C.circuit_val[last_layer_id - 1][j + i * (1 << lg_size)];
        }
    }
}

template<typename FieldT>
FieldT fft_GKR<FieldT>::V_output(const FieldT* output_raw, int r_0_size, int output_size)
{
    FieldT *output;
    output = new FieldT[output_size];
    for(int i = 0; i < output_size; ++i)
        output[i] = output_raw[i];
    for(int i = 0; i < r_0_size; ++i)
    {
        for(int j = 0; j < (output_size >> 1); ++j)
        {
            output[j] = (output[j << 1] * one_minus_r_0[i] + output[j << 1 | 1] * r_0[i]);
        }
        output_size >>= 1;
    }
    FieldT res = output[0];
    delete[] output;
    return res;
}

template<typename FieldT>
quadratic_poly<FieldT> fft_GKR<FieldT>::sumcheck_phase1_update(FieldT previous_random, int current_bit, int total_uv)
{
    quadratic_poly<FieldT> ret = quadratic_poly<FieldT>(FieldT(0), FieldT(0), FieldT(0));
    for(int i = 0; i < (total_uv >> 1); ++i)
    {
        FieldT zero_value, one_value;
        int g_zero = i << 1, g_one = i << 1 | 1;
        if(current_bit == 0)
        {
            V_mult_add[i].b = V_mult_add[g_zero].b;
            V_mult_add[i].a = V_mult_add[g_one].b - V_mult_add[i].b;

            addV_array[i].b = addV_array[g_zero].b;
            addV_array[i].a = addV_array[g_one].b - addV_array[i].b;

            add_mult_sum[i].b = add_mult_sum[g_zero].b;
            add_mult_sum[i].a = add_mult_sum[g_one].b - add_mult_sum[i].b;

        }
        else
        {
            V_mult_add[i].b = (V_mult_add[g_zero].a * previous_random + V_mult_add[g_zero].b);
            V_mult_add[i].a = (V_mult_add[g_one].a * previous_random + V_mult_add[g_one].b - V_mult_add[i].b);

            addV_array[i].b = (addV_array[g_zero].a * previous_random + addV_array[g_zero].b);
            addV_array[i].a = (addV_array[g_one].a * previous_random + addV_array[g_one].b - addV_array[i].b);

            add_mult_sum[i].b = (add_mult_sum[g_zero].a * previous_random + add_mult_sum[g_zero].b);
            add_mult_sum[i].a = (add_mult_sum[g_one].a * previous_random + add_mult_sum[g_one].b - add_mult_sum[i].b);

        }
        ret.a = (ret.a + add_mult_sum[i].a * V_mult_add[i].a);
        ret.b = (ret.b + add_mult_sum[i].a * V_mult_add[i].b + add_mult_sum[i].b * V_mult_add[i].a
                 + addV_array[i].a);
        ret.c = (ret.c + add_mult_sum[i].b * V_mult_add[i].b
                 + addV_array[i].b);
    }
    return ret;
}

template<typename FieldT>
quadratic_poly<FieldT> fft_GKR<FieldT>::sumcheck_phase2_update(FieldT previous_random, int current_bit, int total_uv)
{
    quadratic_poly<FieldT> ret = quadratic_poly<FieldT>(FieldT(0), FieldT(0), FieldT(0));
    for(int i = 0; i < (total_uv >> 1); ++i)
    {
        int g_zero = i << 1, g_one = i << 1 | 1;
        if(current_bit == 0)
        {
            V_mult_add[i].b = V_mult_add[g_zero].b;
            V_mult_add[i].a = V_mult_add[g_one].b - V_mult_add[i].b;

            addV_array[i].b = addV_array[g_zero].b;
            addV_array[i].a = addV_array[g_one].b - addV_array[i].b;

            add_mult_sum[i].b = add_mult_sum[g_zero].b;
            add_mult_sum[i].a = add_mult_sum[g_one].b - add_mult_sum[i].b;
        }
        else
        {

            V_mult_add[i].b = (V_mult_add[g_zero].a * previous_random + V_mult_add[g_zero].b);
            V_mult_add[i].a = (V_mult_add[g_one].a * previous_random + V_mult_add[g_one].b - V_mult_add[i].b);

            addV_array[i].b = (addV_array[g_zero].a * previous_random + addV_array[g_zero].b);
            addV_array[i].a = (addV_array[g_one].a * previous_random + addV_array[g_one].b - addV_array[i].b);

            add_mult_sum[i].b = (add_mult_sum[g_zero].a * previous_random + add_mult_sum[g_zero].b);
            add_mult_sum[i].a = (add_mult_sum[g_one].a * previous_random + add_mult_sum[g_one].b - add_mult_sum[i].b);
        }

        ret.a = (ret.a + add_mult_sum[i].a * V_mult_add[i].a);
        ret.b = (ret.b + add_mult_sum[i].a * V_mult_add[i].b
                 +	add_mult_sum[i].b * V_mult_add[i].a
                 + addV_array[i].a);
        ret.c = (ret.c + add_mult_sum[i].b * V_mult_add[i].b
                 + addV_array[i].b);
    }
    return ret;
}

template<typename FieldT>
bool fft_GKR<FieldT>::addition_layer(FieldT *c_val, int num_poly, FieldT &alpha_beta_sum)
{
    /*
     * addition layer init phase 1
     */
    FieldT zero = FieldT(0);
    beta_g_r0_fhalf[0] = alpha_list[0];
    beta_g_r1_fhalf[0] = beta_list[0];
    beta_g_r0_shalf[0] = FieldT(1);
    beta_g_r1_shalf[0] = FieldT(1);
    int length_g = libff::log2(64);
    int total_uv = poly_size * num_poly;
    int first_half = length_g >> 1, second_half = length_g - first_half;

    for(int i = 0; i < first_half; ++i)
    {
        for(int j = 0; j < (1 << i); ++j)
        {
            beta_g_r0_fhalf[j | (1 << i)] = beta_g_r0_fhalf[j] * r_0[i];
            beta_g_r0_fhalf[j] = beta_g_r0_fhalf[j] * one_minus_r_0[i];
            beta_g_r1_fhalf[j | (1 << i)] = beta_g_r1_fhalf[j] * r_1[i];
            beta_g_r1_fhalf[j] = beta_g_r1_fhalf[j] * one_minus_r_1[i];
        }
    }

    for(int i = 0; i < second_half; ++i)
    {
        for(int j = 0; j < (1 << i); ++j)
        {
            beta_g_r0_shalf[j | (1 << i)] = beta_g_r0_shalf[j] * r_0[i + first_half];
            beta_g_r0_shalf[j] = beta_g_r0_shalf[j] * one_minus_r_0[i + first_half];
            beta_g_r1_shalf[j | (1 << i)] = beta_g_r1_shalf[j] * r_1[i + first_half];
            beta_g_r1_shalf[j] = beta_g_r1_shalf[j] * one_minus_r_1[i + first_half];
        }
    }
    for(int i = 0; i < total_uv; ++i)
    {
        V_mult_add[i] = c_val[i];

        addV_array[i].a = zero;
        addV_array[i].b = zero;
        add_mult_sum[i].a = zero;
        add_mult_sum[i].b = zero;
    }

    int mask_fhalf = (1 << first_half) - 1;
    for(int i = 0; i < num_poly; ++i)
    {
        auto tmp = (beta_g_r0_fhalf[i & mask_fhalf] * beta_g_r0_shalf[i >> first_half]
                    + beta_g_r1_fhalf[i & mask_fhalf] * beta_g_r1_shalf[i >> first_half]);
        for(int j = i * poly_size; j < (i + 1) * poly_size; ++j)
        {
            add_mult_sum[j].b = (add_mult_sum[j].b + tmp);
        }
    }

    //update phase 1
    int log_uv = libff::log2(poly_size) + libff::log2(num_poly);
    FieldT previous_random = FieldT(0);

    for(int i = 0; i < log_uv; ++i)
    {
        quadratic_poly<FieldT> poly = sumcheck_phase1_update(previous_random, i, (1 << (log_uv - i)));
        addition_layer_poly.push_back(poly);
        proof_size += sizeof(quadratic_poly<FieldT>);
        previous_random = r_u_list[0][i];
        if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
        {
            return 0;
        }
        alpha_beta_sum = poly.eval(r_u_list[0][i]);
    }

    v_u_list[0] = V_mult_add[0].eval(previous_random);
    first_half = log_uv >> 1, second_half = log_uv - first_half;

    // update r_0,r_1
    for(int i = 0; i < log_uv; ++i)
    {
        r_0[i] = r_u_list[0][i];
        r_1[i] = r_v_list[0][i];
        one_minus_r_0[i] = FieldT::one() - r_0[i];
        one_minus_r_1[i] = FieldT::one() - r_1[i];
    }

    alpha_beta_sum = alpha_list[0] * v_u_list[0];
    return 1;
}

template<typename FieldT>
bool fft_GKR<FieldT>::addition_layer_verify(FieldT &alpha_beta_sum)
{
    //std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    //int max_bit_length = lg_size + 6;
    //int half_length = (max_bit_length >> 1) + 1;
    //std::vector<FieldT> beta_g_r0_fhalf, beta_g_r0_shalf, beta_g_r1_fhalf, beta_g_r1_shalf, beta_u_fhalf, beta_u_shalf;
    //beta_g_r0_fhalf.resize((lg_size+6) >> 1)+1);

    FieldT zero = FieldT(0);
    beta_g_r0_fhalf[0] = alpha_list[0];
    beta_g_r1_fhalf[0] = beta_list[0];
    beta_g_r0_shalf[0] = FieldT(1);
    beta_g_r1_shalf[0] = FieldT(1);

    int length_g = libff::log2(64);
    int total_uv = poly_size * num_poly;
    int first_half = length_g >> 1, second_half = length_g - first_half;
    int log_uv = libff::log2(poly_size) + libff::log2(num_poly);


    for(int i = 0; i < first_half; ++i)
    {
        for(int j = 0; j < (1 << i); ++j)
        {
            beta_g_r0_fhalf[j | (1 << i)] = beta_g_r0_fhalf[j] * r_0[i];
            beta_g_r0_fhalf[j] = beta_g_r0_fhalf[j] * one_minus_r_0[i];
            beta_g_r1_fhalf[j | (1 << i)] = beta_g_r1_fhalf[j] * r_1[i];
            beta_g_r1_fhalf[j] = beta_g_r1_fhalf[j] * one_minus_r_1[i];
        }
    }

    for(int i = 0; i < second_half; ++i)
    {
        for(int j = 0; j < (1 << i); ++j)
        {
            beta_g_r0_shalf[j | (1 << i)] = beta_g_r0_shalf[j] * r_0[i + first_half];
            beta_g_r0_shalf[j] = beta_g_r0_shalf[j] * one_minus_r_0[i + first_half];
            beta_g_r1_shalf[j | (1 << i)] = beta_g_r1_shalf[j] * r_1[i + first_half];
            beta_g_r1_shalf[j] = beta_g_r1_shalf[j] * one_minus_r_1[i + first_half];
        }
    }

    //update phase 1 : check the sumcheck
    FieldT previous_random = FieldT(0);
    for(int i = 0; i < log_uv; ++i)
    {
        previous_random = r_u_list[0][i];
        if(addition_layer_poly[i].eval(0) + addition_layer_poly[i].eval(1) != alpha_beta_sum)
        {
            return 0;
        }
        alpha_beta_sum = addition_layer_poly[i].eval(r_u_list[0][i]);
    }

    int first_g_half = (length_g >> 1);
    int mask_g_fhalf = (1 << (length_g >> 1)) - 1;
    FieldT summation_val = FieldT(0);
    int log_g = libff::log2(num_poly);
    for(int i = 0; i < num_poly; ++i)
    {
        auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf] * beta_g_r0_shalf[i >> first_g_half]
                      + beta_g_r1_fhalf[i & mask_g_fhalf] * beta_g_r1_shalf[i >> first_g_half]);
        auto tmp_u = FieldT(1);
        for(int j = 0; j < log_g; ++j)
        {
            if((i >> j) & 1)
            {
                tmp_u = tmp_u * r_u_list[0][log_uv - log_g + j];
            }
            else
                tmp_u = tmp_u * (FieldT(1) - r_u_list[0][log_uv - log_g + j]);
        }
        summation_val = summation_val + tmp_g * tmp_u;
    }
    //std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    //v_time += time_span.count();

    for(int i = 0; i < log_uv; ++i)
    {
        r_0[i] = r_u_list[0][i];
        r_1[i] = r_v_list[0][i];
    }
    if(alpha_beta_sum != summation_val * v_u_list[0])
        return 0;
    alpha_beta_sum = alpha_list[0] * v_u_list[0];
    return 1;
}



template<typename FieldT>
bool fft_GKR<FieldT>::mult_layer(FieldT *c_val, int num_poly, FieldT &alpha_beta_sum)
{

    /*
     * mult layer init phase 1
     */
    FieldT zero = FieldT(0);
    beta_g_r0_fhalf[0] = alpha_list[0];
    beta_g_r1_fhalf[0] = beta_list[0];
    beta_g_r0_shalf[0] = FieldT(1);
    beta_g_r1_shalf[0] = FieldT(1);
    int length_g = libff::log2(poly_size * num_poly);
    int total_uv = poly_size;
    int first_half = length_g >> 1, second_half = length_g - first_half;

    for(int i = 0; i < first_half; ++i)
    {
        for(int j = 0; j < (1 << i); ++j)
        {
            beta_g_r0_fhalf[j | (1 << i)] = beta_g_r0_fhalf[j] * r_0[i];
            beta_g_r0_fhalf[j] = beta_g_r0_fhalf[j] * one_minus_r_0[i];
            beta_g_r1_fhalf[j | (1 << i)] = beta_g_r1_fhalf[j] * r_1[i];
            beta_g_r1_fhalf[j] = beta_g_r1_fhalf[j] * one_minus_r_1[i];
        }
    }

    for(int i = 0; i < second_half; ++i)
    {
        for(int j = 0; j < (1 << i); ++j)
        {
            beta_g_r0_shalf[j | (1 << i)] = beta_g_r0_shalf[j] * r_0[i + first_half];
            beta_g_r0_shalf[j] = beta_g_r0_shalf[j] * one_minus_r_0[i + first_half];
            beta_g_r1_shalf[j | (1 << i)] = beta_g_r1_shalf[j] * r_1[i + first_half];
            beta_g_r1_shalf[j] = beta_g_r1_shalf[j] * one_minus_r_1[i + first_half];
        }
    }
    for(int i = 0; i < total_uv; ++i)
    {
        V_mult_add[i] = c_val[i];

        addV_array[i].a = zero;
        addV_array[i].b = zero;
        add_mult_sum[i].a = zero;
        add_mult_sum[i].b = zero;
    }

    int mask_fhalf = (1 << first_half) - 1;
    FieldT *x_eval = new FieldT[num_poly];
    for(int i = 0; i < num_poly; ++i)
        x_eval[i] = FieldT(1);
    for(int i = 0; i < poly_size; ++i)
    {
        for(int j = 0; j < num_poly; ++j)
        {
            int g_id = j * poly_size + i;
            auto tmp = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                        + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);
            add_mult_sum[i].b = (add_mult_sum[i].b + tmp * x_eval[j]);
            x_eval[j] = x_eval[j] * eval_points[j];
        }
    }

    //update phase 1
    int log_uv = libff::log2(poly_size);
    FieldT previous_random = FieldT(0);

    for(int i = 0; i < log_uv; ++i)
    {
        // poly 代表要发送给verifier的sumcheck得到的二次多项式
        quadratic_poly<FieldT> poly = sumcheck_phase1_update(previous_random, i, (1 << (log_uv - i)));
        mult_layer_poly.push_back(poly);
        proof_size += sizeof(quadratic_poly<FieldT>);
        previous_random = r_u_list[1][i];
        if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
        {
            return 0;
        }
        alpha_beta_sum = poly.eval(r_u_list[1][i]);
    }
    v_u_list[1] = V_mult_add[0].eval(previous_random);

    //update r_0,r_1
    for(int i = 0; i < log_uv; ++i)
    {
        r_0[i] = r_u_list[1][i];
        r_1[i] = r_v_list[1][i];
        one_minus_r_0[i] = FieldT::one() - r_0[i];
        one_minus_r_1[i] = FieldT::one() - r_1[i];
    }

    alpha_beta_sum = alpha_list[0] * v_u_list[1];
    return 1;
}

template<typename FieldT>
bool fft_GKR<FieldT>::mult_layer_verify(FieldT &alpha_beta_sum)
{
    //update phase 1
    int log_uv = libff::log2(poly_size);
    FieldT previous_random = FieldT(0);

    for(int i = 0; i < log_uv; ++i)
    {
        // poly 代表要发送给verifier的sumcheck得到的二次多项式
        previous_random = r_u_list[1][i];
        if(mult_layer_poly[i].eval(0) + mult_layer_poly[i].eval(1) != alpha_beta_sum)
        {
            return 0;
        }
        alpha_beta_sum = mult_layer_poly[i].eval(r_u_list[1][i]);
    }

    //std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    auto summation_mult = FieldT(0);
    int length_g = libff::log2(poly_size * num_poly);
    int log_g = length_g;
    for(int i = 0; i < num_poly; ++i)
    {
        auto tmp_g_0 = alpha_list[0];
        auto tmp_g_1 = beta_list[0];
        int lg_num_poly = libff::log2(num_poly);
        for(int j = 0; j < lg_num_poly; ++j)
        {
            if((i >> j) & 1)
            {
                tmp_g_0 = tmp_g_0 * r_0[log_g - lg_num_poly + j];
                tmp_g_1 = tmp_g_1 * r_1[log_g - lg_num_poly + j];
            }
            else
            {
                tmp_g_0 = tmp_g_0 * (FieldT(1) - r_0[log_g - lg_num_poly + j]);
                tmp_g_1 = tmp_g_1 * (FieldT(1) - r_1[log_g - lg_num_poly + j]);
            }
        }
        auto tmp_u_0 = FieldT(1), tmp_u_1 = FieldT(1);
        auto x = eval_points[i];
        for(int j = 0; j < log_uv; ++j)
        {
            tmp_u_0 = tmp_u_0 * (r_0[j] * r_u_list[1][j] * x + (FieldT(1) - r_0[j]) * (FieldT(1) - r_u_list[1][j]));
            tmp_u_1 = tmp_u_1 * (r_1[j] * r_u_list[1][j] * x + (FieldT(1) - r_1[j]) * (FieldT(1) - r_u_list[1][j]));
            x = x * x;
        }
        summation_mult = summation_mult + tmp_g_0 * tmp_u_0 + tmp_g_1 * tmp_u_1;
    }
    //std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    //v_time += time_span.count();
    for(int i = 0; i < log_uv; ++i)
    {
        r_0[i] = r_u_list[1][i];
        r_1[i] = r_v_list[1][i];
    }
    delete[] eval_points;
    eval_points = NULL;
    if(alpha_beta_sum != summation_mult * v_u_list[1])
        return 0;
    alpha_beta_sum = alpha_list[0] * v_u_list[1];
    return 1;
}

template<typename FieldT>
bool fft_GKR<FieldT>::intermediate_layer(FieldT &alpha_beta_sum)
{
    if(invert)
        alpha_beta_sum = alpha_beta_sum * FieldT(poly_size);
    assert(FieldT(poly_size) * inv_n == FieldT(1));
    return true;
}

template<typename FieldT>
bool fft_GKR<FieldT>::ifft_gkr(circuit<FieldT> &C, FieldT &alpha_beta_sum)
{
    std::vector<FieldT> rot_mul;
    rot_mul.resize(62,0);
    rot_mul[0] = invert? inv_rou : rou;
    for(int i = 1; i < 62; ++i)
    {
        rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1];
    }
    size_t starting_depth = 0;
    for(size_t dep = 0, index = 0; dep < poly_dimension; ++dep, ++index)
    {
        int blk_size = 1 << (lg_size - dep);
        int half_blk_size = blk_size >> 1;
        int cur = starting_depth + (poly_dimension - dep); // from poly_dimension to 1
        int pre = cur - 1;
        {
            //printf("cur = %lu, dep = %lu\n", cur, dep);
            FieldT zero = FieldT(0);
            beta_g_r0_fhalf[0] = alpha_list[index]; //mark dep =
            beta_g_r1_fhalf[0] = beta_list[index]; //
            beta_g_r0_shalf[0] = FieldT(1);
            beta_g_r1_shalf[0] = FieldT(1);
            int length_g = libff::log2(poly_size);
            int total_uv = poly_size;
            int first_half = length_g >> 1, second_half = length_g - first_half;

            for(int i = 0; i < first_half; ++i)
            {
                for(int j = 0; j < (1 << i); ++j)
                {
                    beta_g_r0_fhalf[j | (1 << i)] = beta_g_r0_fhalf[j] * r_0[i];
                    beta_g_r0_fhalf[j] = beta_g_r0_fhalf[j] * one_minus_r_0[i];
                    beta_g_r1_fhalf[j | (1 << i)] = beta_g_r1_fhalf[j] * r_1[i];
                    beta_g_r1_fhalf[j] = beta_g_r1_fhalf[j] * one_minus_r_1[i];
                }
            }

            for(int i = 0; i < second_half; ++i)
            {
                for(int j = 0; j < (1 << i); ++j)
                {
                    beta_g_r0_shalf[j | (1 << i)] = beta_g_r0_shalf[j] * r_0[i + first_half];
                    beta_g_r0_shalf[j] = beta_g_r0_shalf[j] * one_minus_r_0[i + first_half];
                    beta_g_r1_shalf[j | (1 << i)] = beta_g_r1_shalf[j] * r_1[i + first_half];
                    beta_g_r1_shalf[j] = beta_g_r1_shalf[j] * one_minus_r_1[i + first_half];
                }
            }
            for(int i = 0; i < total_uv; ++i)
            {
                V_mult_add[i] = C.circuit_val[pre][i];

                addV_array[i].a = zero;
                addV_array[i].b = zero;
                add_mult_sum[i].a = zero;
                add_mult_sum[i].b = zero;
            }

            int mask_fhalf = (1 << first_half) - 1;
            FieldT x = FieldT(1);
            for(int k = 0; k < blk_size / 2; ++k)
            {
                for(int j = 0; j < (1 << dep); ++j)
                {
                    int g_id = k << (dep) | j;
                    int u_id = k << (dep + 1) | j;
                    int v_id = k << (dep + 1) | (1 << dep) | j;
                    auto tmp = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                                + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);

                    add_mult_sum[u_id].b = add_mult_sum[u_id].b + tmp;
                    addV_array[u_id].b = addV_array[u_id].b + tmp * x * C.circuit_val[pre][v_id];


                    g_id = (k + blk_size / 2) << (dep) | j;
                    u_id = k << (dep + 1) | j;
                    v_id = k << (dep + 1) | (1 << dep) | j;
                    tmp = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                           + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);
                    add_mult_sum[u_id].b = add_mult_sum[u_id].b + tmp;
                    addV_array[u_id].b = addV_array[u_id].b - tmp * x * C.circuit_val[pre][v_id];
                }
                x = x * rot_mul[dep];
            }

            //update phase 1
            int log_uv = libff::log2(poly_size);
            FieldT previous_random = FieldT(0);
            std::vector<quadratic_poly<FieldT>> temp1;
            for(int i = 0; i < log_uv; ++i)
            {
                quadratic_poly<FieldT> poly = sumcheck_phase1_update(previous_random, i, (1 << (log_uv - i)));
                temp1.push_back(poly);
                proof_size += sizeof(quadratic_poly<FieldT>);
                previous_random = r_u_list[index+2][i];
                if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
                {
                    printf("phase1 error! cur = %lu  i = %d", cur, i);
                    assert(false);
                    return 0;
                }
                alpha_beta_sum = poly.eval(r_u_list[index+2][i]);
            }
            fft_layer_phase1_poly.push_back(temp1);
            v_u_list[index+2] = V_mult_add[0].eval(previous_random);

            //phase 2 init
            beta_u_fhalf[0] = FieldT(1);
            beta_u_shalf[0] = FieldT(1);
            for(int i = 0; i < first_half; ++i)
            {
                for(int j = 0; j < (1 << i); ++j)
                {
                    beta_u_fhalf[j | (1 << i)] = beta_u_fhalf[j] * r_u_list[index+2][i];
                    beta_u_fhalf[j] = beta_u_fhalf[j] * (FieldT::one() - r_u_list[index+2][i]);
                }
            }

            for(int i = 0; i < second_half; ++i)
            {
                for(int j = 0; j < (1 << i); ++j)
                {
                    beta_u_shalf[j | (1 << i)] = beta_u_shalf[j] * r_u_list[index+2][i + first_half];
                    beta_u_shalf[j] = beta_u_shalf[j] * (FieldT::one() - r_u_list[index+2][i + first_half]);
                }
            }
            for(int i = 0; i < total_uv; ++i)
            {
                V_mult_add[i] = C.circuit_val[pre][i];

                addV_array[i].a = zero;
                addV_array[i].b = zero;
                add_mult_sum[i].a = zero;
                add_mult_sum[i].b = zero;
            }

            x = FieldT(1);
            for(int k = 0; k < blk_size / 2; ++k)
            {
                for(int j = 0; j < (1 << dep); ++j)
                {
                    int g_id = k << (dep) | j;
                    int u_id = k << (dep + 1) | j;
                    int v_id = k << (dep + 1) | (1 << dep) | j;
                    auto tmp_g = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                                  + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);
                    auto tmp_u = (beta_u_fhalf[u_id & mask_fhalf] * beta_u_shalf[u_id >> first_half]);
                    auto tmp_g_u = tmp_g * tmp_u;
                    add_mult_sum[v_id].b = add_mult_sum[v_id].b + tmp_g_u * x;
                    addV_array[v_id].b = addV_array[v_id].b + tmp_g_u * v_u_list[index+2];


                    g_id = (k + blk_size / 2) << (dep) | j;
                    u_id = k << (dep + 1) | j;
                    v_id = k << (dep + 1) | (1 << dep) | j;
                    tmp_g = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                             + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);
                    tmp_u = (beta_u_fhalf[u_id & mask_fhalf] * beta_u_shalf[u_id >> first_half]);
                    tmp_g_u = tmp_g * tmp_u;
                    add_mult_sum[v_id].b = add_mult_sum[v_id].b - tmp_g_u * x;
                    addV_array[v_id].b = addV_array[v_id].b + tmp_g_u * v_u_list[index+2];
                }
                x = x * rot_mul[dep];
            }


            //phase 2 update, I think phase 2 can be removed
            previous_random = FieldT(0);
            std::vector<quadratic_poly<FieldT>> temp2;
            for(int i = 0; i < log_uv; ++i)
            {
                quadratic_poly<FieldT> poly = sumcheck_phase2_update(previous_random, i, (1 << (log_uv - i)));
                temp2.push_back(poly);
                proof_size += sizeof(quadratic_poly<FieldT>);
                previous_random = r_v_list[index+2][i];
                if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
                {
                    printf("phase2 error! cur = %lu  i = %d",cur,i);
                    assert(false);
                    return 0;
                }
                alpha_beta_sum = poly.eval(r_v_list[index+2][i]);
            }
            fft_layer_phase2_poly.push_back(temp2);
            v_v_list[index+2] = V_mult_add[0].eval(previous_random);

            //update r_0,r_1
            for(int i = 0; i < log_uv; ++i)
            {
                r_0[i] = r_u_list[index+2][i];
                r_1[i] = r_v_list[index+2][i];
                one_minus_r_0[i] = FieldT::one() - r_u_list[index+2][i];
                one_minus_r_1[i] = FieldT::one() - r_v_list[index+2][i];
            }
            alpha_beta_sum = alpha_list[index+1] * v_u_list[index+2] + beta_list[index+1] * v_v_list[index+2];
        }
    }
    return 1;
}

template<typename FieldT>
bool fft_GKR<FieldT>::ifft_gkr_verify(FieldT &alpha_beta_sum)
{
    std::vector<FieldT> rot_mul;
    rot_mul.resize(62,0);
    rot_mul[0] = invert? inv_rou : rou;
    for(int i = 1; i < 62; ++i)
    {
        rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1];
    }
    size_t starting_depth = 0;
    int log_uv = libff::log2(poly_size);
    std::vector<quadratic_poly<FieldT>> *temp1,*temp2;
    for(size_t dep = 0, index = 0; dep < poly_dimension; ++dep, ++index)
    {
        int blk_size = 1 << (lg_size - dep);
        int half_blk_size = blk_size >> 1;
        int cur = starting_depth + (poly_dimension - dep);
        int pre = cur - 1;
        {
            //update phase 1
            //assert(alpha_beta_sum == (alpha_list[0] * v_u_list[1]));
            FieldT previous_random = FieldT(0);
            temp1 = &fft_layer_phase1_poly[index];
            for(int i = 0; i < log_uv; ++i)
            {
                quadratic_poly<FieldT> poly = temp1->at(i);
                previous_random = r_u_list[index+2][i];
                if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
                {
                    printf("phase1 error! cur = %lu  i = %d", cur, i);
                    assert(false);
                    return 0;
                }
                alpha_beta_sum = poly.eval(r_u_list[index+2][i]);
            }

            //phase 2 update, I think phase 2 can be removed
            previous_random = FieldT(0);
            temp2 = &fft_layer_phase2_poly[index];
            for(int i = 0; i < log_uv; ++i)
            {
                quadratic_poly<FieldT> poly = temp2->at(i);
                previous_random = r_v_list[index+2][i];
                if(poly.eval(0) + poly.eval(1) != alpha_beta_sum)
                {
                    printf("phase2 error! cur = %lu  i = %d",cur,i);
                    assert(false);
                    return 0;
                }
                alpha_beta_sum = poly.eval(r_v_list[index+2][i]);
            }

            //verifier's computation
            //printf("cur = %lu, dep = %lu\n", cur, dep);
            //std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
            auto summation_beta_u_0_A = FieldT(1), summation_beta_v_0_A = FieldT(1);
            auto summation_beta_u_1_A = FieldT(1), summation_beta_v_1_A = FieldT(1);
            auto x = rot_mul[dep];
            size_t log_k = libff::log2(blk_size / 2);
            size_t log_j = dep;
            summation_beta_u_0_A = (FieldT(1) - r_0[log_uv - 1]) * (FieldT(1) - r_u_list[index+2][log_j]) * r_v_list[index+2][log_j] * alpha_list[index];  //mark
            summation_beta_u_1_A = (FieldT(1) - r_1[log_uv - 1]) * (FieldT(1) - r_u_list[index+2][log_j]) * r_v_list[index+2][log_j] * beta_list[index];
            summation_beta_v_0_A = (FieldT(1) - r_0[log_uv - 1]) * (FieldT(1) - r_u_list[index+2][log_j]) * r_v_list[index+2][log_j] * alpha_list[index];
            summation_beta_v_1_A = (FieldT(1) - r_1[log_uv - 1]) * (FieldT(1) - r_u_list[index+2][log_j]) * r_v_list[index+2][log_j] * beta_list[index];

            auto summation_beta_u_0_B = r_0[log_uv - 1] * (FieldT(1) - r_u_list[index+2][log_j]) * r_v_list[index+2][log_j] * alpha_list[index]; //mark
            auto summation_beta_u_1_B = r_1[log_uv - 1] * (FieldT(1) - r_u_list[index+2][log_j]) * r_v_list[index+2][log_j] * beta_list[index];
            auto summation_beta_v_0_B = r_0[log_uv - 1] * (FieldT(1) - r_u_list[index+2][log_j]) * r_v_list[index+2][log_j] * alpha_list[index];
            auto summation_beta_v_1_B = r_1[log_uv - 1] * (FieldT(1) - r_u_list[index+2][log_j]) * r_v_list[index+2][log_j] * beta_list[index];
            for(int i = 0; i < log_k; ++i)
            {
                summation_beta_u_0_A = summation_beta_u_0_A * (r_0[log_j + i] * r_u_list[index+2][log_j + 1 + i] * r_v_list[index+2][log_j + 1 + i] + (FieldT(1) - r_0[log_j + i]) * (FieldT(1) - r_u_list[index+2][log_j + 1 + i]) * (FieldT(1) - r_v_list[index+2][log_j + 1 + i]));
                summation_beta_u_1_A = summation_beta_u_1_A * (r_1[log_j + i] * r_u_list[index+2][log_j + 1 + i] * r_v_list[index+2][log_j + 1 + i] + (FieldT(1) - r_1[log_j + i]) * (FieldT(1) - r_u_list[index+2][log_j + 1 + i]) * (FieldT(1) - r_v_list[index+2][log_j + 1 + i]));
                summation_beta_v_0_A = summation_beta_v_0_A * (r_0[log_j + i] * r_u_list[index+2][log_j + 1 + i] * r_v_list[index+2][log_j + 1 + i] * x + (FieldT(1) - r_0[log_j + i]) * (FieldT(1) - r_u_list[index+2][log_j + 1 + i]) * (FieldT(1) - r_v_list[index+2][log_j + 1 + i]));
                summation_beta_v_1_A = summation_beta_v_1_A * (r_1[log_j + i] * r_u_list[index+2][log_j + 1 + i] * r_v_list[index+2][log_j + 1 + i] * x + (FieldT(1) - r_1[log_j + i]) * (FieldT(1) - r_u_list[index+2][log_j + 1 + i]) * (FieldT(1) - r_v_list[index+2][log_j + 1 + i]));

                summation_beta_u_0_B = summation_beta_u_0_B * (r_0[log_j + i] * r_u_list[index+2][log_j + 1 + i] * r_v_list[index+2][log_j + 1 + i] + (FieldT(1) - r_0[log_j + i]) * (FieldT(1) - r_u_list[index+2][log_j + 1 + i]) * (FieldT(1) - r_v_list[index+2][log_j + 1 + i]));
                summation_beta_u_1_B = summation_beta_u_1_B * (r_1[log_j + i] * r_u_list[index+2][log_j + 1 + i] * r_v_list[index+2][log_j + 1 + i] + (FieldT(1) - r_1[log_j + i]) * (FieldT(1) - r_u_list[index+2][log_j + 1 + i]) * (FieldT(1) - r_v_list[index+2][log_j + 1 + i]));
                summation_beta_v_0_B = summation_beta_v_0_B * ((FieldT(1) - r_0[log_j + i]) * (FieldT(1) - r_u_list[index+2][log_j + 1 + i]) * (FieldT(1) - r_v_list[index+2][log_j + 1 + i]) + r_0[log_j + i] * r_u_list[index+2][log_j + 1 + i] * r_v_list[index+2][log_j + 1 + i] * x);
                summation_beta_v_1_B = summation_beta_v_1_B * ((FieldT(1) - r_1[log_j + i]) * (FieldT(1) - r_u_list[index+2][log_j + 1 + i]) * (FieldT(1) - r_v_list[index+2][log_j + 1 + i]) + r_1[log_j + i] * r_u_list[index+2][log_j + 1 + i] * r_v_list[index+2][log_j + 1 + i] * x);
                x = x * x;
            }
            for(int i = 0; i < log_j; ++i)
            {
                summation_beta_u_0_A = summation_beta_u_0_A * (r_0[i] * r_u_list[index+2][i] * r_v_list[index+2][i] + (FieldT(1) - r_0[i]) * (FieldT(1) - r_u_list[index+2][i]) * (FieldT(1) - r_v_list[index+2][i]));
                summation_beta_u_1_A = summation_beta_u_1_A * (r_1[i] * r_u_list[index+2][i] * r_v_list[index+2][i] + (FieldT(1) - r_1[i]) * (FieldT(1) - r_u_list[index+2][i]) * (FieldT(1) - r_v_list[index+2][i]));
                summation_beta_v_0_A = summation_beta_v_0_A * (r_0[i] * r_u_list[index+2][i] * r_v_list[index+2][i] + (FieldT(1) - r_0[i]) * (FieldT(1) - r_u_list[index+2][i]) * (FieldT(1) - r_v_list[index+2][i]));
                summation_beta_v_1_A = summation_beta_v_1_A * (r_1[i] * r_u_list[index+2][i] * r_v_list[index+2][i] + (FieldT(1) - r_1[i]) * (FieldT(1) - r_u_list[index+2][i]) * (FieldT(1) - r_v_list[index+2][i]));

                summation_beta_u_0_B = summation_beta_u_0_B * (r_0[i] * r_u_list[index+2][i] * r_v_list[index+2][i] + (FieldT(1) - r_0[i]) * (FieldT(1) - r_u_list[index+2][i]) * (FieldT(1) - r_v_list[index+2][i]));
                summation_beta_u_1_B = summation_beta_u_1_B * (r_1[i] * r_u_list[index+2][i] * r_v_list[index+2][i] + (FieldT(1) - r_1[i]) * (FieldT(1) - r_u_list[index+2][i]) * (FieldT(1) - r_v_list[index+2][i]));
                summation_beta_v_0_B = summation_beta_v_0_B * (r_0[i] * r_u_list[index+2][i] * r_v_list[index+2][i] + (FieldT(1) - r_0[i]) * (FieldT(1) - r_u_list[index+2][i]) * (FieldT(1) - r_v_list[index+2][i]));
                summation_beta_v_1_B = summation_beta_v_1_B * (r_1[i] * r_u_list[index+2][i] * r_v_list[index+2][i] + (FieldT(1) - r_1[i]) * (FieldT(1) - r_u_list[index+2][i]) * (FieldT(1) - r_v_list[index+2][i]));
            }

            if(alpha_beta_sum != (summation_beta_u_0_A + summation_beta_u_1_A + summation_beta_u_0_B + summation_beta_u_1_B) * v_u_list[index+2] + (summation_beta_v_0_A + summation_beta_v_1_A - summation_beta_v_0_B - summation_beta_v_1_B) * v_v_list[index+2])
                return 0;

            //std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
            //std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
            //v_time += time_span.count();

            for(int i = 0; i < log_uv; ++i)
            {
                r_0[i] = r_u_list[index+2][i];
                r_1[i] = r_v_list[index+2][i];
            }
        }
        alpha_beta_sum = alpha_list[index+1] * v_u_list[index+2] + beta_list[index+1] * v_v_list[index+2];
    }
    return 1;
}

template<typename FieldT>
bool fft_GKR<FieldT>::extension_gkr()
{
    size_t lg_size_ = libff::log2(poly_size);
    for(int i = 1; i <= lg_size_; ++i)
    {
        for(int j = 0; j < i; ++j)
            proof_size += sizeof(quadratic_poly<FieldT>);
    }
    return true;
}

template<typename FieldT>
int fft_GKR<FieldT>::engage_gkr(const std::vector<FieldT> &poly_coeffs, std::vector<FieldT> &result)
{
    circuit<FieldT> C;
    build_circuit(C,poly_coeffs,result);

    /* init array */
    init_array(lg_size + 6);
    r_0.resize(lg_size + 10);
    one_minus_r_0.resize(lg_size + 10);
    r_1.resize(lg_size + 10);
    one_minus_r_1.resize(lg_size + 10);

    refresh_randomness(r_0, one_minus_r_0, lg_size + 10);
    refresh_randomness(r_1, one_minus_r_1, lg_size + 10);

    std::vector<FieldT> r0,r1,one_minus_r0(one_minus_r_0),one_minus_r1(one_minus_r_1);
    r0.assign(r_0.begin(), r_0.end());
    r1.assign(r_1.begin(), r_1.end());

    alpha_beta_sum0 =  V_output(C.circuit_val[C.circuit_val.size() - 1], libff::log2(C.size[C.size.size() - 1]), C.size[C.size.size() - 1]);
    FieldT alpha_beta_sum = alpha_beta_sum0;
    bool res = true;

    //poly eval part
    //last layer, only addition
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    res &= addition_layer(C.circuit_val[C.circuit_val.size() - 2],  64, alpha_beta_sum);
    res &= mult_layer(C.circuit_val[C.circuit_val.size() - 3], 64, alpha_beta_sum);

    //printf("inter_layer_id: %lu\n", C.circuit_val.size() - 4);
    res &= intermediate_layer(alpha_beta_sum);
    assert(res);
    //ifft part
    res &= ifft_gkr(C,alpha_beta_sum);
    assert(res);
    //extension part
    res &= extension_gkr();
    assert(res);
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    p_time += time_span.count();

    r_0.assign(r0.begin(),r0.end());
    r_1.assign(r1.begin(),r1.end());
    one_minus_r_0.assign(one_minus_r0.begin(), one_minus_r0.end());
    one_minus_r_1.assign(one_minus_r1.begin(), one_minus_r1.end());

    destruction_array();
    return res;
}

template<typename FieldT>
void fft_GKR<FieldT>::init_array(int max_bit_length, bool prover)
{
    if(prover) {
        add_mult_sum = new linear_poly<FieldT>[1 << max_bit_length];
        V_mult_add = new linear_poly<FieldT>[1 << max_bit_length];
        addV_array = new linear_poly<FieldT>[1 << max_bit_length];

        // init alpha,beta list

        random_field_vector<FieldT>(alpha_list,poly_dimension + 1);
        random_field_vector<FieldT>(beta_list,poly_dimension + 1);
        alpha_list[0] = FieldT::one();
        beta_list[0] = FieldT::zero();

        // init r_u, r_v
        for(int i = 0; i < poly_dimension + 2; i++)
        {
            std::vector<FieldT> temp1, temp2;
            random_field_vector<FieldT>(temp1,lg_size + 10);
            random_field_vector<FieldT>(temp2,lg_size + 10);
            r_u_list.push_back(temp1);
            r_v_list.push_back(temp2);
        }

        // init v_v, v_u
        v_v_list.resize(poly_dimension + 2);
        v_u_list.resize(poly_dimension + 2);

    }
    int half_length = (max_bit_length >> 1) + 1;
    beta_g_r0_fhalf = new FieldT[(1 << half_length)];
    beta_g_r0_shalf = new FieldT[(1 << half_length)];
    beta_g_r1_fhalf = new FieldT[(1 << half_length)];
    beta_g_r1_shalf = new FieldT[(1 << half_length)];
    beta_u_fhalf = new FieldT[(1 << half_length)];
    beta_u_shalf = new FieldT[(1 << half_length)];
}

template<typename FieldT>
bool fft_GKR<FieldT>::destruction_array(bool prover)
{
    if(prover){
        delete[] add_mult_sum; add_mult_sum = NULL;
        delete[] V_mult_add;  V_mult_add = NULL;
        delete[] addV_array; addV_array = NULL;
    }
    delete[] beta_g_r0_fhalf; beta_g_r0_fhalf = NULL;
    delete[] beta_g_r0_shalf; beta_g_r0_shalf = NULL;
    delete[] beta_g_r1_fhalf; beta_g_r1_fhalf = NULL;
    delete[] beta_g_r1_shalf; beta_g_r1_shalf = NULL;
    delete[] beta_u_fhalf; beta_u_fhalf = NULL;
    delete[] beta_u_shalf; beta_u_shalf = NULL;
}


template<typename FieldT>
bool fft_GKR<FieldT>::prover_compute(const std::vector<FieldT> &poly_coeffs, std::vector<FieldT> &result)
{
    bool res = engage_gkr(poly_coeffs,result);
    assert(res);
    proof_size += sizeof(FieldT) * (poly_dimension + 2) * 2; // the size of v_u,v_v
    //printf("prover compute successfully. size of poly = %lu\n", poly_size);

    return res;
}

template<typename FieldT>
bool fft_GKR<FieldT>::verifier_predicate()
{
    bool res = true;
    init_array(lg_size + 6,0);
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    FieldT alpha_beta_sum = alpha_beta_sum0;
    /* add layer */
    res &= addition_layer_verify(alpha_beta_sum);
    assert(res);
    /* mult layer */
    res &= mult_layer_verify(alpha_beta_sum);

    res &= intermediate_layer(alpha_beta_sum);

    assert(res);
    res &= ifft_gkr_verify(alpha_beta_sum);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    v_time = time_span.count();

    if(!res)
        fprintf(stderr, "Error, fft gkr failed\n");
    //printf("vtime %f\nproof size %d\nptime %f\n", v_time, proof_size, p_time);

    destruction_array(false);
    return res;
}

template<typename FieldT>
void fft_GKR<FieldT>::refresh_randomness(FieldT *r, FieldT *one_minus_r, int size)
{
    for(int i = 0; i < size; ++i)
    {
        r[i] = FieldT::random_element();
        one_minus_r[i] = FieldT(1) - r[i];
    }
}

template<typename FieldT>
void fft_GKR<FieldT>::refresh_randomness(std::vector<FieldT> &r, std::vector<FieldT> &one_minus_r, int size)
{
    for(int i = 0; i < size; ++i)
    {
        r[i] = FieldT::random_element();
        one_minus_r[i] = FieldT::one() - r[i];
    }
}

template<typename FieldT>
void random_field_vector(std::vector<FieldT> &v, size_t n)
{
    assert(v.size() <= n);
    if(v.size() < n)
        v.resize(n, FieldT::zero());
    for(size_t i = 0; i < n; i++)
    {
        v[i] = FieldT::random_element();
    }
}

}


#endif //LIGERO_FFT_CIRCUIT_GKR_TCC
