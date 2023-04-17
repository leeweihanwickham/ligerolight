/** 本部分的主要作用在于计算编码后矩阵的行长、列长，计算交互阶段和访问阶段的次数等 **/
#include <python2.7/Python.h>
#include "ligero_iop.hpp"

namespace ligero {
template<typename FieldT>
ligero_iop_parameters<FieldT>::ligero_iop_parameters(const std::size_t security_parameter,
                                                     const std::size_t RS_extra_dimensions,
                                                     const std::size_t RS_col_extra_dimensions,
                                                     const std::vector<std::size_t> localization_parameter_array,
                                                     // 描述了待编码前的谕示矩阵的横纵比
                                                     const float height_width_ratio,
                                                     const bool make_zk,
                                                     const field_subset_type domain_type,
                                                     // 乘法门的约束数目，也就是x,y,z的长度
                                                     const std::size_t num_constraints,
                                                     // 总变量数，也就是w的长度，P_x w = x, P_y w = y, P_z w = z
                                                     const std::size_t num_variables,
                                                     const std::size_t parallel_opt) :
    security_parameter_(security_parameter),
    RS_extra_dimensions_(RS_extra_dimensions),
    RS_col_extra_dimensions_(RS_col_extra_dimensions),
    localization_parameter_array_(localization_parameter_array),
    height_width_ratio_(height_width_ratio),
    make_zk_(make_zk),
    domain_type_(domain_type),
    num_constraints_(num_constraints),
    num_variables_(num_variables),
    parallel_opt_(parallel_opt)
{
    this->set_soundness_parameters();
    this->print();
}

template<typename FieldT>
void ligero_iop_parameters<FieldT>::set_soundness_parameters()
{

    /** Ligero++的可靠性误差计算方式如下，当 e=d/3, d=|L|-k+1，如果k=|H|-1，则d=|L|-|H|+2
     * 共有4次 Testing Interleaved RS code，4d/|F| + 4*(1-e/|L|)^{t}
     * 3次 Testing Linear Constraints，1/|F|+((e+k+l)/|L|)^{t}，这里l是编码前的向量长度，l \leq k-1
     * 1次 Testing Quadratic Constraints, 1/|F|+((e+2k)/|L|)^{t}
     *
     * Besides, there will be IPA soundness error
     * IPA interactive: 1/|F|^{sigma2}
     * IPA query: security_parameter/(IPA_RS_Extra_dim)
     *
     * 因此，总可靠性误差可认为是
     * epsilon = (4d+4))/(|F|^{sigma})+ 4*(1-e/|L|)^{t}+4*((2k+e)/|L|)^{t}
     * 其中，k=|H|-1, sigma表示interactive phase 的重复次数
     * 对于严格零知识性的版本，度会增大，|H|->|H'|，|H'| \geq |H|+t，我的算法是确定矩阵行宽比时先确定有效行长再加上t，从而找到|H'|
     *
     * ligero++的可靠性误差共来自于3部分，因此每部分的可靠性应该达到 2^{-security_parameter - 2}
     * **/

    const size_t interactive_soundness_bits = this->security_parameter_ + 1;
    const size_t query_soundness_bits = this->security_parameter_ + 1;
    /** Parameterize direct LDT and RS-encoded Ligero's queries. */
    this->configure_encoded_ligero_params(this->num_variables_, this->num_constraints_);
    this->set_queries(query_soundness_bits);
    this->set_encoded_ligero_interactions(interactive_soundness_bits);
}

template<typename FieldT>
void ligero_iop_parameters<FieldT>::set_encoded_ligero_interactions(const size_t interactive_soundness_bits)
{
    /** Encoded protocol's interactive soundness error is (4d + 3 + 1 / F)^(num_interactive_repetitions)
     *  Consequently the bits of security is:
     *      -Interactive_soundness_bits = (num_interactive_repetitions) * log(4d+4 / |F|)
     *  Rearranging, we get:
     *      num_interactive_repetitions = ceil(-Interactive_soundness_bits / log(4d+4 / |F|))
     *  We compute the logarithm using the identity:
     *      log(1 / |F|) = -log(|F|)
     */

    // 其实相当于表示了log_2 |F|
    const long double field_size_bits = (long double)(libff::soundness_log_of_field_size_helper<FieldT>(FieldT::zero()));
    const long double code_distance = ((1ull << this->codeword_domain_dim_) - (1ull << this->systematic_domain_dim_) + 1) / 2;
    const long double denominator = log2( 4 * code_distance + 4 )-field_size_bits;
    const long double num_interactive_repetitions = ceil(-1.0 * interactive_soundness_bits / denominator);
    this->encoded_ligero_params_.num_interaction_phase_repetitions_ =
        std::max<size_t>(1, size_t(num_interactive_repetitions));

}

template<typename FieldT>
void ligero_iop_parameters<FieldT>::set_queries(const size_t query_soundness_bits)
{
    /** queries in ligerolight is all columns*/
    this->encoded_ligero_params_.num_query_phase_repetitions_ = (1ull << this->codeword_domain_dim_);
}


template<typename FieldT>
void ligero_iop_parameters<FieldT>::configure_encoded_ligero_params(const size_t num_variables,
                                                                    const size_t num_constraints)
{

    std::size_t num_vars = this->num_variables_ + 1;
    // 表示待编码前的谕示矩阵行长，ceil表示取整，根据长宽比直接算出的行长
    // TODO zk
    std::size_t systematic_domain_size = (std::size_t) ceil(sqrt(num_vars / this->height_width_ratio_)) ;
    // 由于要运行FFT，行长要是2的某次幂；此外，还要增加t
    this->systematic_domain_dim_ = libff::log2(systematic_domain_size);
    // TODO zk
    this->codeword_domain_dim_ = this->systematic_domain_dim_ + this->RS_extra_dimensions_;

    systematic_domain_size = 1ull << this->systematic_domain_dim_;
    //表示U_w谕示矩阵的高度，每行最多放|H|-t个有效元素
    // TODO zk
    this->num_oracles_input_ =  ceil((float) num_vars) / (systematic_domain_size) + 1;
    //this->num_oracles_input_ = (size_t) ceil(((float) num_vars) / (systematic_domain_size));
    // 表示P_x，P_y, P_z等矩阵的宽度
    std::size_t matrix_width = (systematic_domain_size) * this->num_oracles_input_;

    //num_oracle_vectors表示U_x,U_y,U_z等谕示矩阵的高度，类似的，也最多只能放|H|-t个有效元素
    std::size_t matrix_height = num_constraints;

    if (matrix_height % (systematic_domain_size ) != 0)
    {
        matrix_height += (systematic_domain_size)
                - matrix_height % (systematic_domain_size);
    }
    this->num_oracle_vectors_ = matrix_height / (systematic_domain_size);
    // TODO zk

    std::size_t col_systematic_domain_size =
            libff::round_to_next_power_of_2(this->num_oracle_vectors_ + this->num_oracles_input_);
    std::size_t col_systematic_domain_dim = libff::log2(col_systematic_domain_size);
    this->col_systematic_domain_dim_ = col_systematic_domain_dim;

    std::size_t total_localizations = 0;
    for (size_t i = 0; i < localization_parameter_array_.size(); i++) {
        total_localizations += localization_parameter_array_[i];
    }
    assert(total_localizations <= (col_systematic_domain_dim+1));

    const std::size_t summation_degree_bound = 2 * (1ull << col_systematic_domain_dim) - 1;
    const long double field_size_bits = (long double)(libff::soundness_log_of_field_size_helper<FieldT>(FieldT::zero()));
    std::size_t IPA_inter_repetition =
            ceil(( this->security_parameter_ - 1 + libff::log2(summation_degree_bound) ) / field_size_bits);
    const std::size_t IPA_num_interaction_repetitions = std::max<size_t> (1, IPA_inter_repetition);

    const std::size_t IPA_num_query_repetitions = ceil((this->security_parameter_ / this->RS_col_extra_dimensions_)) + 1;

    this->encoded_ligero_params_.make_zk_ = this->make_zk_;
    this->encoded_ligero_params_.field_subset_type_ = this->domain_type_;
    this->encoded_ligero_params_.matrix_width_ = matrix_width;
    this->encoded_ligero_params_.matrix_height_ = matrix_height;
    this->encoded_ligero_params_.num_oracles_input_ = this->num_oracles_input_;
    this->encoded_ligero_params_.num_oracles_vectors_ = this->num_oracle_vectors_;
    this->encoded_ligero_params_.col_systematic_domain_dim_ = col_systematic_domain_dim;
    this->encoded_ligero_params_.localization_parameter_array = this->localization_parameter_array_;
    this->encoded_ligero_params_.IPA_num_interaction_repetitions_ = IPA_num_interaction_repetitions;
    this->encoded_ligero_params_.IPA_num_query_repetitions_ = IPA_num_query_repetitions;
}

template<typename FieldT>
size_t ligero_iop_parameters<FieldT>::systematic_domain_dim() const
{
    return this->systematic_domain_dim_;
}
template<typename FieldT>
size_t ligero_iop_parameters<FieldT>::RS_extra_dimensions() const
{
    return this->RS_extra_dimensions_;
}

template<typename FieldT>
std::size_t ligero_iop_parameters<FieldT>::col_RS_extra_dimensions() const {
    return this->RS_col_extra_dimensions_;
}

template<typename FieldT>
bool ligero_iop_parameters<FieldT>::make_zk() const
{
    return this->make_zk_;
}
template<typename FieldT>
field_subset_type ligero_iop_parameters<FieldT>::domain_type() const
{
    return this->domain_type_;
}

template<typename FieldT>
void ligero_iop_parameters<FieldT>::print() const
{
    printf("\nLigero IOP parameters\n");
    libff::print_indent(); printf("* target security parameter = %zu\n", this->security_parameter_);

    libff::print_indent(); printf("* encoded ligero interactions = %lu\n",
        this->encoded_ligero_params_.num_interaction_phase_repetitions_);
    libff::print_indent(); printf("* encoded ligero queries = %lu\n",
        this->encoded_ligero_params_.num_query_phase_repetitions_);
    libff::print_indent(); printf("* RS extra dimensions = %zu\n", this->RS_extra_dimensions_);
    libff::print_indent(); printf("* IPA interactions = %lu\n",
                                  this->encoded_ligero_params_.IPA_num_interaction_repetitions_);
    libff::print_indent(); printf("* IPA queries = %lu\n",
                                  this->encoded_ligero_params_.IPA_num_query_repetitions_);
    libff::print_indent(); printf("* col_systematic_domain_dim = %zu\n", this->col_systematic_domain_dim_);
    libff::print_indent(); printf("* IPA extra dimensions = %zu\n", this->RS_col_extra_dimensions_);
    libff::print_indent(); printf("* localization_parameter_array = ");
    for (std::size_t i = 0; i < this->localization_parameter_array_.size(); i ++)
    {
        libff::print_indent(); std::cout << " " << this->localization_parameter_array_[i] ;
    }
    std::cout << std::endl;
    libff::print_indent(); printf("* systematic domain dim = %zu\n", this->systematic_domain_dim_);
    libff::print_indent(); printf("* codeword domain dim = %zu\n", this->codeword_domain_dim_);
    libff::print_indent(); printf("* num oracles for input = %zu\n", this->num_oracles_input_);
    libff::print_indent(); printf("* num oracle vectors = %zu\n", this->num_oracle_vectors_);
    libff::print_indent(); printf("* make zk = %s\n", (this->make_zk_ ? "true" : "false"));
    libff::print_indent(); printf("* domain type = %s\n", field_subset_type_names[this->domain_type_]);
}


template<typename FieldT>
ligero_iop<FieldT>::ligero_iop(iop_protocol<FieldT> &IOP,
                               const r1cs_constraint_system<FieldT> &constraint_system,
                               const ligero_iop_parameters<FieldT> &parameters) :
    IOP_(IOP),
    constraint_system_(constraint_system),
    parameters_(parameters)
{
    const size_t systematic_domain_size = 1ull << parameters.systematic_domain_dim();
    const size_t codeword_domain_dim = parameters.systematic_domain_dim() + this->parameters_.RS_extra_dimensions();
    const size_t codeword_domain_size = 1ull << codeword_domain_dim;

    const std::size_t extended_systematic_domain_size = systematic_domain_size << 1;
    FieldT shift;
    if (parameters_.domain_type() == multiplicative_coset_type)
    {
        shift = FieldT::multiplicative_generator;
    }
    else
    {
        shift = FieldT(codeword_domain_size);
    }
    assert(this->parameters_.col_systematic_domain_dim_ == this->parameters_.encoded_ligero_params_.col_systematic_domain_dim_);
    const size_t col_systematic_domain_dim = this->parameters_.encoded_ligero_params_.col_systematic_domain_dim_;
    const std::size_t col_systematic_domain_size = 1ull << col_systematic_domain_dim;
    const size_t col_codeword_domain_dim = col_systematic_domain_dim + this->parameters_.col_RS_extra_dimensions();
    const size_t col_codeword_domain_size= 1ull<<col_codeword_domain_dim;

    this->codeword_domain_ = field_subset<FieldT>(codeword_domain_size);
    this->RS_col_codeword_domain_=field_subset<FieldT>(col_codeword_domain_size,FieldT(col_codeword_domain_size));

    field_subset<FieldT> systematic_domain(systematic_domain_size);
    field_subset<FieldT> extended_systematic_domain(extended_systematic_domain_size);
    field_subset<FieldT> RS_col_systematic_domain(col_systematic_domain_size);

    domain_handle codeword_domain_handle = this->IOP_.register_domain(this->codeword_domain_);
    domain_handle systematic_domain_handle = this->IOP_.register_domain(systematic_domain);
    domain_handle RS_col_systematic_domain_handle = this->IOP_.register_domain(RS_col_systematic_domain);
    domain_handle RS_col_codeword_domain_handle = this->IOP_.register_domain(this->RS_col_codeword_domain_);
    domain_handle extended_systematic_domain_handle = this->IOP_.register_domain(extended_systematic_domain);

    this->protocol_ = std::make_shared<interleaved_r1cs_protocol<FieldT> >(this->IOP_,
                                                                           codeword_domain_handle,
                                                                           systematic_domain_handle,
                                                                           extended_systematic_domain_handle,
                                                                           RS_col_systematic_domain_handle,
                                                                           RS_col_codeword_domain_handle,
                                                                           constraint_system,
                                                                           this->parameters_.encoded_ligero_params_,
                                                                           this->parameters_.parallel_opt_,
                                                                           this->parameters_.RS_col_extra_dimensions_);
}

template<typename FieldT>
void ligero_iop<FieldT>::register_interactions(){}

template<typename FieldT>
void ligero_iop<FieldT>::register_queries()
{
    // 不能删 这个是生成打开的列
    this->protocol_->register_queries();
}

template<typename FieldT>
float ligero_iop<FieldT>::produce_oracle(const r1cs_primary_input<FieldT> &primary_input,
                                                const r1cs_auxiliary_input<FieldT> &auxiliary_input,
                                                const std::map<std::size_t,FieldT> public_input) {

    if (this->parameters_.make_zk()) {
        this->protocol_->submit_blinding_vector_oracles();
    }

    float mt_time = this->protocol_->submit_witness_oracles(primary_input, auxiliary_input, public_input);
    return mt_time;
}

template<typename FieldT>
float ligero_iop<FieldT>::produce_proof(const std::map<std::size_t,FieldT> public_input) {

    float time = this->protocol_->calculate_and_submit_proof(public_input);

    return time;
}

template<typename FieldT>
bool ligero_iop<FieldT>::verifier_predicate(const std::map<std::size_t,FieldT> public_input)
{
    libff::enter_block("Check Interleaved R1CS verifier predicate");
    bool decision = this->protocol_->verifier_predicate(public_input);
    libff::leave_block("Check Interleaved R1CS verifier predicate");
    return decision;
}

template<typename FieldT>
void ligero_iop<FieldT>::submit_random_blinding_vector(const oracle_handle_ptr &handle)
{
    polynomial<FieldT> random_poly = polynomial<FieldT>::random_polynomial(this->systematic_domain_size_);
    std::vector<FieldT> random_vector = FFT_over_field_subset<FieldT>(random_poly.coefficients(), this->codeword_domain_);
    oracle<FieldT> random_vector_oracle(random_vector);
    this->IOP_.submit_oracle(handle, std::move(random_vector_oracle));
}
}// ligero
