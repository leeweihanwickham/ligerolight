#include <cmath>
#include <stdexcept>
#include <sys/time.h>
#include <libff/common/profiling.hpp>
#include "ligero/algebra/fft.hpp"
#include <libff/common/utils.hpp>
#include "ligero/protocols/ldt/fri/fri_ldt.hpp"

namespace ligero {

template<typename FieldT>
interleaved_r1cs_protocol<FieldT>::interleaved_r1cs_protocol(
    iop_protocol<FieldT> &IOP,
    const domain_handle &codeword_domain_handle,
    const domain_handle &systematic_domain_handle,
    const domain_handle &extended_systematic_domain_handle,
    const domain_handle &RS_col_systematic_domain_handle,
    const domain_handle &RS_col_codeword_domain_handle,
    const r1cs_constraint_system<FieldT> &constraint_system,
    const encoded_ligero_parameters &parameters,
    const std::size_t parallel_opt,
    const std::size_t IPA_RS_extra_dimensions) :
    IOP_(IOP),
    constraint_system_(constraint_system),
    parameters_(parameters),
    codeword_domain_handle_(codeword_domain_handle),
    systematic_domain_handle_(systematic_domain_handle),
    extended_systematic_domain_handle_(extended_systematic_domain_handle),
    RS_col_systematic_domain_handle_(RS_col_systematic_domain_handle),
    RS_col_codeword_domain_handle_(RS_col_codeword_domain_handle),
    parallel_opt_(parallel_opt),
    IPA_RS_extra_dimensions_(IPA_RS_extra_dimensions)
{
    this->num_queries_ = parameters.num_query_phase_repetitions_;

    this->num_interactions_ = parameters.num_interaction_phase_repetitions_;
    this->make_zk_ = parameters.make_zk_;
    this->field_subset_type_ = parameters.field_subset_type_;

    this->matrix_width_ = parameters.matrix_width_;
    this->matrix_height_ = parameters.matrix_height_;
    this->num_oracles_input_ = parameters.num_oracles_input_;
    this->num_oracles_vectors_ = parameters.num_oracles_vectors_;
//    this->MT_security_parameter = parameters.MT_security_parameter;

    this->codeword_domain_ = this->IOP_.get_domain(this->codeword_domain_handle_);
    this->systematic_domain_ = this->IOP_.get_domain(this->systematic_domain_handle_);
    this->extended_systematic_domain_ = this->IOP_.get_domain(this->extended_systematic_domain_handle_);
    this->RS_col_codeword_domain_=this->IOP_.get_domain(this->RS_col_codeword_domain_handle_);
    this->RS_col_systematic_domain_=this->IOP_.get_domain(this->RS_col_systematic_domain_handle_);

    this->codeword_domain_size_ = this->codeword_domain_.num_elements();
    this->systematic_domain_size_ = this->systematic_domain_.num_elements();
    this->extended_systematic_domain_size_ = this->extended_systematic_domain_.num_elements();
    this->RS_col_codeword_domain_size_=this->RS_col_codeword_domain_.num_elements();
    this->RS_col_systematic_domain_size_=this->RS_col_systematic_domain_.num_elements();

    this->encoding_independence_ = 3;

    /* Get R1CS matrices */
    this->A_matrix_ = this->constraint_system_.A_matrix();
    this->B_matrix_ = this->constraint_system_.B_matrix();
    this->C_matrix_ = this->constraint_system_.C_matrix();
    this->value_matrix_for_commitment.resize(this->num_oracles_input_ + 3 * this->num_oracles_vectors_);
    /* Add extra rows */
    this->A_matrix_.resize(this->matrix_height_);
    this->B_matrix_.resize(this->matrix_height_);
    this->C_matrix_.resize(this->matrix_height_);

    this->lincheck_A_.reset(new interleaved_lincheck_ot_protocol<FieldT>(this->IOP_,
                                                                         this->codeword_domain_handle_,
                                                                         this->systematic_domain_handle_,
                                                                         this->extended_systematic_domain_handle_,
                                                                         this->RS_col_systematic_domain_handle_,
                                                                         this->num_oracles_input_,//m2
                                                                         this->num_oracles_vectors_,//m1
                                                                         this->num_queries_,
                                                                         this->num_interactions_,
                                                                         this->make_zk_,
                                                                         this->field_subset_type_,
                                                                         this->A_matrix_,
                                                                         this->parallel_opt_));
    this->lincheck_B_.reset(new interleaved_lincheck_ot_protocol<FieldT>(this->IOP_,
                                                                         this->codeword_domain_handle_,
                                                                         this->systematic_domain_handle_,
                                                                         this->extended_systematic_domain_handle_,
                                                                         this->RS_col_systematic_domain_handle_,
                                                                         this->num_oracles_input_,
                                                                         this->num_oracles_vectors_,
                                                                         this->num_queries_,
                                                                         this->num_interactions_,
                                                                         this->make_zk_,
                                                                         this->field_subset_type_,
                                                                         this->B_matrix_,
                                                                         this->parallel_opt_));
    this->lincheck_C_.reset(new interleaved_lincheck_ot_protocol<FieldT>(this->IOP_,
                                                                         this->codeword_domain_handle_,
                                                                         this->systematic_domain_handle_,
                                                                         this->extended_systematic_domain_handle_,
                                                                         this->RS_col_systematic_domain_handle_,
                                                                         this->num_oracles_input_,
                                                                         this->num_oracles_vectors_,
                                                                         this->num_queries_,
                                                                         this->num_interactions_,
                                                                         this->make_zk_,
                                                                         this->field_subset_type_,
                                                                         this->C_matrix_,
                                                                         this->parallel_opt_));

    this->IRScheck_w_.reset(new interleaved_RScode<FieldT>(this->IOP_,
                            this->codeword_domain_handle_,
                            this->systematic_domain_handle_,
                            this->RS_col_systematic_domain_handle_,
                            this->num_oracles_input_,
                            this->num_queries_,
                            this->num_interactions_,
                            this->make_zk_,
                            this->field_subset_type_));

    this->IRScheck_x_.reset(new interleaved_RScode<FieldT>(this->IOP_,
                            this->codeword_domain_handle_,
                            this->systematic_domain_handle_,
                            this->RS_col_systematic_domain_handle_,
                            this->num_oracles_vectors_,
                            this->num_queries_,
                            this->num_interactions_,
                            this->make_zk_,
                            this->field_subset_type_));

    this->IRScheck_y_.reset(new interleaved_RScode<FieldT>(this->IOP_,
                            this->codeword_domain_handle_,
                            this->systematic_domain_handle_,
                            this->RS_col_systematic_domain_handle_,
                            this->num_oracles_vectors_,
                            this->num_queries_,
                            this->num_interactions_,
                            this->make_zk_,
                            this->field_subset_type_));

    this->IRScheck_z_.reset(new interleaved_RScode<FieldT>(this->IOP_,
                            this->codeword_domain_handle_,
                            this->systematic_domain_handle_,
                            this->RS_col_systematic_domain_handle_,
                            this->num_oracles_vectors_,
                            this->num_queries_,
                            this->num_interactions_,
                            this->make_zk_,
                            this->field_subset_type_));


    this->rowcheck_.reset(new interleaved_rowcheck_protocol<FieldT>(this->IOP_,
                                                                    this->codeword_domain_handle_,//L
                                                                    this->systematic_domain_handle_,//H
                                                                    this->extended_systematic_domain_handle_,//2H
                                                                    this->RS_col_systematic_domain_handle_,
                                                                    this->num_oracles_input_,
                                                                    this->num_oracles_vectors_,//m1
                                                                    this->num_queries_,
                                                                    this->num_interactions_,
                                                                    this->make_zk_,
                                                                    this->field_subset_type_));

    if (this->make_zk_)
    {
//        this->lincheck_blinding_vector.resize(this->num_interactions_);
        this->rowcheck_blinding_vector.resize(this->num_interactions_);
    }
}


template<typename FieldT>
void interleaved_r1cs_protocol<FieldT>::register_queries()
{
    std::vector<random_query_position_handle> query_position_handles;
    query_position_handles.resize(this->num_queries_);

    //TODO change
    for (std::size_t i = 0; i < this->num_queries_; ++i){
        this->query_set.emplace_back(i);
    }

//    for (std::size_t i = 0; i < this->num_queries_; ++i)
//    {
//        bool is_repeat = true;
//        std::size_t val;
//        while(is_repeat){
//            query_position_handles[i] = this->IOP_.register_random_query_position(this->codeword_domain_handle_);
//            val = this->IOP_.obtain_query_position(query_position_handles[i]);
//            std::vector<std::size_t>::iterator it;
//            it = find(this->query_set.begin(),this->query_set.end(),val);
//            if(it == this->query_set.end()){
//                is_repeat = false;
//            }
//        }
//        this->query_set.emplace_back(val);
//    }
}


/* Proving */
template<typename FieldT>
float interleaved_r1cs_protocol<FieldT>::submit_witness_oracles(const std::vector<FieldT> &primary_input,
                                                               const std::vector<FieldT> &auxiliary_input,//construct complete witness
                                                               const std::map<std::size_t, FieldT> &public_input)//public witness
{
    libff::enter_block("Submit witness oracles");

    libff::enter_block("Generate extended witness and auxiliary witness");

    /* construct z = (1, v, w) */
    std::vector<FieldT> extended_witness = std::vector<FieldT>(1, FieldT(1));
    extended_witness.insert(extended_witness.end(),
                            primary_input.begin(),
                            primary_input.end());
    extended_witness.insert(extended_witness.end(),
                            auxiliary_input.begin(),
                            auxiliary_input.end());

    const std::size_t extra_elements = this->matrix_width_ - extended_witness.size();
    for (std::size_t i = 0; i < extra_elements; ++i)
    {
        extended_witness.emplace_back(FieldT(0));
    }

    /* Extended witness, with public input replaced by 0s */
    std::vector<FieldT> auxiliary_only_witness = extended_witness;//w_bar
    typename std::map<std::size_t, FieldT>::const_iterator it;
    for (it = public_input.begin(); it != public_input.end(); it++)
    {
        std::size_t idx = it->first;
        FieldT val = it->second;
        if(auxiliary_only_witness[idx] == val){
            auxiliary_only_witness[idx] = FieldT(0);
        }
	//else public_input has sth wrong?
    }
    libff::leave_block("Generate extended witness and auxiliary witness");

    libff::enter_block("Perform matrix multiplications");
    //attach value to a,b,c
    std::vector<FieldT> a_result_vector;
    for (std::size_t i = 0; i < this->matrix_height_; ++i)
    {
        FieldT sum(0);
        std::map<std::size_t, FieldT> row = this->A_matrix_[i];
        typename std::map<std::size_t, FieldT>::iterator it;
        for (it = row.begin(); it != row.end(); it++)
        {
            std::size_t idx = it->first;
            FieldT val = it->second;
            sum += val * extended_witness[idx];
        }
        a_result_vector.emplace_back(sum);
    }

    std::vector<FieldT> b_result_vector;
    for (std::size_t i = 0; i < this->matrix_height_; ++i)
    {
        FieldT sum(0);
        std::map<std::size_t, FieldT> row = this->B_matrix_[i];
        typename std::map<std::size_t, FieldT>::iterator it;
        for (it = row.begin(); it != row.end(); it++)
        {
            std::size_t idx = it->first;
            FieldT val = it->second;
            sum += val * extended_witness[idx];
        }
        b_result_vector.emplace_back(sum);
    }

    std::vector<FieldT> c_result_vector;
    for (std::size_t i = 0; i < this->matrix_height_; ++i)
    {
        FieldT sum(0);
        std::map<std::size_t, FieldT> row = this->C_matrix_[i];
        typename std::map<std::size_t, FieldT>::iterator it;
        for (it = row.begin(); it != row.end(); it++)
        {
            std::size_t idx = it->first;
            FieldT val = it->second;
            sum += val * extended_witness[idx];
        }
        c_result_vector.emplace_back(sum);
    }
    libff::leave_block("Perform matrix multiplications");
    //submit oracles
    //用来填充的全0向量
    const std::vector<FieldT> PaddingZeroVec(this->codeword_domain_size_,FieldT::zero());
    //将谕示矩阵Uwxy扩展到2^k以编码
    assert(libff::round_to_next_power_of_2(this->num_oracles_vectors_+this->num_oracles_input_)==RS_col_systematic_domain_size_);
    //初始化扩展后的谕示矩阵 此时Uwxyz未编码，是systematic_size 按行存储

    this->Uw_matrix.resize(RS_col_systematic_domain_size_,PaddingZeroVec);
    this->Ux_matrix.resize(RS_col_systematic_domain_size_,PaddingZeroVec);
    this->Uy_matrix.resize(RS_col_systematic_domain_size_,PaddingZeroVec);
    this->Uz_matrix.resize(RS_col_systematic_domain_size_,PaddingZeroVec);

    libff::enter_block("Submit input oracles");
    for (std::size_t i = 0; i < this->num_oracles_input_; ++i)
    {
        // TODO zk
        // arrange w into m_2 * (h-t) matrix W, arrange W to m_2 * h
        // w->U_w
        const std::size_t start = i * (this->systematic_domain_size_ );
        const std::size_t end = start + (this->systematic_domain_size_ );

        std::vector<FieldT> w_row(&auxiliary_only_witness[start], &auxiliary_only_witness[end]);//a row of matrix W

        // TODO zk
        const std::vector<FieldT> w_row_coefficients =
            IFFT_over_field_subset<FieldT>(w_row, this->systematic_domain_);
	    std::vector<FieldT> w_row_vector = FFT_over_field_subset(w_row_coefficients, this->codeword_domain_);
        //TODO 能否只保留matrix_w_commit,去掉this->Uw_matrix？  answer:不能，Uw是按行存储，在正确性检查时用编码前的“编码矩阵”  Uw还没有RS编码，matrix_commit是转置，按列存储，编码后的矩阵
        // 因此这里不计算matrix_w_commit
       //Uw是在填充的0矩阵上面，所以索引从0->num_oracles_input_，而且也不用考虑make_zk向量的位置，因为它也是在Uw下面
        this->Uw_matrix[i]=std::move(w_row_vector);
    }
    assert(this->Uw_matrix[this->num_oracles_input_]==PaddingZeroVec);
    libff::leave_block("Submit input oracles");

    libff::enter_block("Submit vector oracles");
    libff::enter_block("Horizon Encoding");
    for (std::size_t i = 0; i < this->num_oracles_vectors_; ++i)
    {
        // TODO zk
        // arrange a,b,c into m_1 * h matrix A,B,C
        const std::size_t start =  i * (this->systematic_domain_size_);
        const std::size_t end = start + this->systematic_domain_size_ ;

        std::vector<FieldT> a_row(&a_result_vector[start], &a_result_vector[end]);
        std::vector<FieldT> b_row(&b_result_vector[start], &b_result_vector[end]);
        std::vector<FieldT> c_row(&c_result_vector[start], &c_result_vector[end]);
        // TODO zk
        const std::vector<FieldT> a_row_coefficients =
            IFFT_over_field_subset<FieldT>(a_row, this->systematic_domain_);
        std::vector<FieldT> a_row_vector = FFT_over_field_subset(a_row_coefficients, this->codeword_domain_);

        const std::vector<FieldT> b_row_coefficients =
                IFFT_over_field_subset<FieldT>(b_row, this->systematic_domain_);
        std::vector<FieldT> b_row_vector = FFT_over_field_subset(b_row_coefficients, this->codeword_domain_);

        const std::vector<FieldT> c_row_coefficients =
                IFFT_over_field_subset<FieldT>(c_row, this->systematic_domain_);
        std::vector<FieldT> c_row_vector = FFT_over_field_subset(c_row_coefficients, this->codeword_domain_);
//         这里要考虑到Uxyz的位置了
        if(this->make_zk_){
//            给blinding留一行
            this->Ux_matrix[this->num_oracles_input_+1+i]=std::move(a_row_vector);
            this->Uy_matrix[this->num_oracles_input_+1+i]=std::move(b_row_vector);
            this->Uz_matrix[this->num_oracles_input_+1+i]=std::move(c_row_vector);
        } else{
            this->Ux_matrix[this->num_oracles_input_+i]=std::move(a_row_vector);
            this->Uy_matrix[this->num_oracles_input_+i]=std::move(b_row_vector);
            this->Uz_matrix[this->num_oracles_input_+i]=std::move(c_row_vector);
        }
    }

    libff::leave_block("Horizon Encoding");

    libff::enter_block("Vertical Encoding");
    //对blinding矩阵编码
    if (this->make_zk_){
        //编码后按列存储的矩阵

        // do not use class member, use normal variable instead
//        std::vector<std::vector<FieldT>> matrix_Uxyz_Rowcheck_blinding_commit;

        matrix_Uxyz_Rowcheck_blinding_commit.resize(this->codeword_domain_size_);
        //编码前按行存储的矩阵
        // TODO OPT: 不用计算编码前按行存储的矩阵 直接可以给出转置后的矩阵

        rowcheck_blinding_matrix.resize(RS_col_systematic_domain_size_,PaddingZeroVec);
        rowcheck_blinding_matrix[this->num_oracles_input_+this->num_oracles_vectors_+1]=this->rowcheck_blinding_vector[0];

        for(std::size_t k=0;k<this->codeword_domain_size_;k++){

            std::vector<FieldT> rowcheck_blinding_slice(this->RS_col_systematic_domain_size_);
            //矩阵转置
            for(std::size_t j=0;j<this->RS_col_systematic_domain_size_;j++){
                rowcheck_blinding_slice[j]=rowcheck_blinding_matrix[j][k];
            }

            //对列编码
            const std::vector<FieldT> rowcheck_col_coefficients=IFFT_over_field_subset<FieldT>(rowcheck_blinding_slice, this->RS_col_systematic_domain_);
            std::vector<FieldT> rowcheck_col_vector=FFT_over_field_subset(rowcheck_col_coefficients, this->RS_col_codeword_domain_);
            matrix_Uxyz_Rowcheck_blinding_commit[k] = std::move(rowcheck_col_vector);
        }
    }

    // matrix_commit 按列存放 是对Uwxyz的列编码
    // do not use class member, use normal variable instead

//    std::vector<std::vector<FieldT>> matrix_w_commit;
//    std::vector<std::vector<FieldT>> matrix_x_commit;
//    std::vector<std::vector<FieldT>> matrix_y_commit;
//    std::vector<std::vector<FieldT>> matrix_z_commit;

    matrix_w_commit.resize(this->codeword_domain_size_);
    matrix_x_commit.resize(this->codeword_domain_size_);
    matrix_y_commit.resize(this->codeword_domain_size_);
    matrix_z_commit.resize(this->codeword_domain_size_);
    // 对Uwxyz纵向编码
    for(std::size_t i=0;i<this->codeword_domain_size_;i++){
        std::vector<FieldT> w_slice(this->RS_col_systematic_domain_size_);
        std::vector<FieldT> x_slice(this->RS_col_systematic_domain_size_);
        std::vector<FieldT> y_slice(this->RS_col_systematic_domain_size_);
        std::vector<FieldT> z_slice(this->RS_col_systematic_domain_size_);
        for(std::size_t j=0;j<this->RS_col_systematic_domain_size_;j++){
            w_slice[j]=this->Uw_matrix[j][i];
            x_slice[j]=this->Ux_matrix[j][i];
            y_slice[j]=this->Uy_matrix[j][i];
            z_slice[j]=this->Uz_matrix[j][i];
        }
        const std::vector<FieldT> w_col_coefficients=IFFT_over_field_subset<FieldT>(w_slice, this->RS_col_systematic_domain_);
        const std::vector<FieldT> x_col_coefficients=IFFT_over_field_subset<FieldT>(x_slice, this->RS_col_systematic_domain_);
        const std::vector<FieldT> y_col_coefficients=IFFT_over_field_subset<FieldT>(y_slice, this->RS_col_systematic_domain_);
        const std::vector<FieldT> z_col_coefficients=IFFT_over_field_subset<FieldT>(z_slice, this->RS_col_systematic_domain_);
        std::vector<FieldT> w_col_vector=FFT_over_field_subset(w_col_coefficients, this->RS_col_codeword_domain_);
        std::vector<FieldT> x_col_vector=FFT_over_field_subset(x_col_coefficients, this->RS_col_codeword_domain_);
        std::vector<FieldT> y_col_vector=FFT_over_field_subset(y_col_coefficients, this->RS_col_codeword_domain_);
        std::vector<FieldT> z_col_vector=FFT_over_field_subset(z_col_coefficients, this->RS_col_codeword_domain_);
        matrix_w_commit[i]=std::move(w_col_vector);
        matrix_x_commit[i]=std::move(x_col_vector);
        matrix_y_commit[i]=std::move(y_col_vector);
        matrix_z_commit[i]=std::move(z_col_vector);
    }

    libff::leave_block("Vertical Encoding");
    libff::leave_block("Submit vector oracles");
    libff::leave_block("Submit witness oracles");

    libff::enter_block("generate MT root and commitment");

    struct timeval start1,end1;
    gettimeofday(&start1, nullptr);

    // current commit,
    std::vector<std::size_t> merkle_queries=this->merkleTree->get_queries
            (this->parameters_.IPA_num_query_repetitions_,(this->RS_col_codeword_domain_size_>>this->parameters_.localization_parameter_array[0]));

        // construct the matrix to be committed
    std::vector<std::vector<FieldT>> matrix_w_commit_matrix;
    std::vector<std::vector<FieldT>> matrix_x_commit_matrix;
    std::vector<std::vector<FieldT>> matrix_y_commit_matrix;
    std::vector<std::vector<FieldT>> matrix_z_commit_matrix;
    std::vector<std::vector<FieldT>> matrix_blinding_commit_matrix;
    std::vector<std::vector<FieldT>> matrix_together_commit_matrix;

    /**Note: here we support put all the oracle matrices into one matrix
     * This is feasible because queries are all the same**/

    matrix_w_commit_matrix.resize((1ull << this->parameters_.localization_parameter_array[0])* this->codeword_domain_size_);
    matrix_x_commit_matrix.resize((1ull << this->parameters_.localization_parameter_array[0])* this->codeword_domain_size_);
    matrix_y_commit_matrix.resize((1ull << this->parameters_.localization_parameter_array[0])* this->codeword_domain_size_);
    matrix_z_commit_matrix.resize((1ull << this->parameters_.localization_parameter_array[0])* this->codeword_domain_size_);
    matrix_blinding_commit_matrix.resize((1ull << this->parameters_.localization_parameter_array[0])* this->codeword_domain_size_);
    matrix_together_commit_matrix.resize((1ull << this->parameters_.localization_parameter_array[0])* this->codeword_domain_size_ * 5);

    std::size_t set_size = (this->RS_col_codeword_domain_size_ >> this->parameters_.localization_parameter_array[0]);

    for (std::size_t i = 0; i < this->codeword_domain_size_; i ++)
    {
        for ( std::size_t j = 0; j < (1ull << this->parameters_.localization_parameter_array[0]); j ++)
        {
            matrix_w_commit_matrix[i + j*matrix_w_commit.size()].assign(matrix_w_commit[i].begin()+j*set_size,matrix_w_commit[i].begin()+(j+1)*set_size);
            matrix_x_commit_matrix[i + j*matrix_w_commit.size()].assign(matrix_x_commit[i].begin()+j*set_size,matrix_x_commit[i].begin()+(j+1)*set_size);
            matrix_y_commit_matrix[i + j*matrix_w_commit.size()].assign(matrix_y_commit[i].begin()+j*set_size,matrix_y_commit[i].begin()+(j+1)*set_size);
            matrix_z_commit_matrix[i + j*matrix_w_commit.size()].assign(matrix_z_commit[i].begin()+j*set_size,matrix_z_commit[i].begin()+(j+1)*set_size);
            matrix_blinding_commit_matrix[i + j*matrix_w_commit.size()].assign(matrix_Uxyz_Rowcheck_blinding_commit[i].begin()+j*set_size,matrix_Uxyz_Rowcheck_blinding_commit[i].begin()+(j+1)*set_size);
        }
    }

    for (std::size_t i = 0; i < this->codeword_domain_size_ ;i ++)
    {
        matrix_together_commit_matrix[i] = matrix_w_commit_matrix[i];
        matrix_together_commit_matrix[i + this->codeword_domain_size_] = matrix_w_commit_matrix[i + this->codeword_domain_size_];
        matrix_together_commit_matrix[i + this->codeword_domain_size_*2] = matrix_x_commit_matrix[i];
        matrix_together_commit_matrix[i + this->codeword_domain_size_*3] = matrix_x_commit_matrix[i + this->codeword_domain_size_];
        matrix_together_commit_matrix[i + this->codeword_domain_size_*4] = matrix_y_commit_matrix[i];
        matrix_together_commit_matrix[i + this->codeword_domain_size_*5] = matrix_y_commit_matrix[i + this->codeword_domain_size_];
        matrix_together_commit_matrix[i + this->codeword_domain_size_*6] = matrix_z_commit_matrix[i];
        matrix_together_commit_matrix[i + this->codeword_domain_size_*7] = matrix_z_commit_matrix[i + this->codeword_domain_size_];
        matrix_together_commit_matrix[i + this->codeword_domain_size_*8] = matrix_blinding_commit_matrix[i];
        matrix_together_commit_matrix[i + this->codeword_domain_size_*9] = matrix_blinding_commit_matrix[i + this->codeword_domain_size_];
    }


    gettimeofday(&end1, nullptr);
    const float mt_time = (end1.tv_usec-start1.tv_usec)/1000000.0 + end1.tv_sec-start1.tv_sec;

    // current commit
    this->merkleTree.reset(new merkle<FieldT>(
            (this->RS_col_codeword_domain_size_>> this->parameters_.localization_parameter_array[0]),
            merkle_queries,
            true
    ));

    this->lincheck_blinding_commit.resize(1);

//    首先计算承诺路径
    // current commit
    std::vector<std::size_t>auxiliary_pos=this->merkleTree->find_merkle_path_only_index
            (2*(this->RS_col_codeword_domain_size_ >> this->parameters_.localization_parameter_array[0])-1);

    libff::enter_block("zk MT");

    libff::leave_block("zk MT");
    libff::enter_block("oracle MT");
// 对Uwxyz承诺只需一次，不要交互次数

    this->lincheck_blinding_commit[0] = std::move(this->merkleTree->create_merklePar_of_matrix(matrix_together_commit_matrix));

    std::vector<std::vector<FieldT>>().swap(matrix_w_commit_matrix);
    std::vector<std::vector<FieldT>>().swap(matrix_x_commit_matrix);
    std::vector<std::vector<FieldT>>().swap(matrix_y_commit_matrix);
    std::vector<std::vector<FieldT>>().swap(matrix_z_commit_matrix);
    std::vector<std::vector<FieldT>>().swap(matrix_blinding_commit_matrix);
    std::vector<std::vector<FieldT>>().swap(matrix_together_commit_matrix);

    libff::leave_block("oracle MT");
    libff::leave_block("generate MT root and commitment");

    return mt_time;

}

template<typename FieldT>
std::vector<FieldT> interleaved_r1cs_protocol<FieldT>::submit_zero_blinding_vector()
{
    std::vector<FieldT> elems(this->extended_systematic_domain_size_, FieldT(0));
    for (std::size_t i = this->systematic_domain_size_; i < this->systematic_domain_size_ + this->encoding_independence_; ++i)
    {
        const std::size_t index = this->extended_systematic_domain_.reindex_by_subset(
            this->systematic_domain_.dimension(), i);
        elems[index]=FieldT::random_element();
    }

    const std::vector<FieldT> coeffs = IFFT_over_field_subset<FieldT>(elems, this->extended_systematic_domain_);
    const std::vector<FieldT> vector = FFT_over_field_subset<FieldT>(coeffs, this->codeword_domain_);
    return vector;
}


template<typename FieldT>
void interleaved_r1cs_protocol<FieldT>::submit_blinding_vector_oracles()
{
    assert(this->make_zk_);
//    this->lincheck_blinding_vector[0]=std::move(this->submit_zero_sum_blinding_vector());
    this->rowcheck_blinding_vector[0]=std::move(this->submit_zero_blinding_vector());
    for (std::size_t i = 1; i < this->num_interactions_; ++i)
    {
//        this->lincheck_blinding_vector[i]=this->lincheck_blinding_vector[0];
        this->rowcheck_blinding_vector[i]=this->rowcheck_blinding_vector[0];
    }
}


template<typename FieldT>
float interleaved_r1cs_protocol<FieldT>::calculate_and_submit_proof(const std::map<std::size_t,FieldT> &public_input)
{

    struct timeval start2,end2;
    gettimeofday(&start2, nullptr);

    this->lincheck_ABC_randomcombinations.resize(this->num_interactions_);
    for(std::size_t h=0;h<this->num_interactions_;h++){
        this->lincheck_ABC_randomcombinations[h]= std::move(random_FieldT_vector<FieldT>(this->A_matrix_.size()));
    }

    this->IRS_random_combinations.resize(this->num_interactions_);
    for(std::size_t h=0;h<this->num_interactions_;h++){
        this->IRS_random_combinations[h]=std::move(random_FieldT_vector<FieldT>(this->Uw_matrix.size()));
    }

    gettimeofday(&end2, nullptr);
    const float random_time = (end2.tv_usec-start2.tv_usec)/1000000.0 + end2.tv_sec-start2.tv_sec ;

    libff::enter_block("Calculating and submitting proof");
    /* construct additional input as (1,v,0) */
    std::vector<FieldT> additional_input = std::vector<FieldT>(this->matrix_width_, FieldT(0));//w_bar
    typename std::map<std::size_t, FieldT>::const_iterator it;
    for (it = public_input.begin(); it != public_input.end(); it++)
    {
        std::size_t idx = it->first;
        FieldT val = it->second;
        additional_input[idx] = val;
	//else public_input has sth wrong?
    }//v_bar

    libff::enter_block("Calculating encoding matrix of additional_input");
    const std::vector<FieldT> PaddingZeroVec(this->codeword_domain_size_,FieldT::zero());
    this->additional_matrix.resize(this->RS_col_systematic_domain_size_,PaddingZeroVec);
    for (std::size_t i = 0; i < this->num_oracles_input_; ++i)
    {
        const std::size_t vecstart = i * (this->systematic_domain_size_ );
        const std::size_t vecend = vecstart + (this->systematic_domain_size_ );

        std::vector<FieldT> row(&additional_input[vecstart], &additional_input[vecend]);

        const std::vector<FieldT> row_coefficients =
                IFFT_over_field_subset<FieldT>(row, this->systematic_domain_);
        std::vector<FieldT> row_vector =
                FFT_over_field_subset<FieldT>(row_coefficients, this->codeword_domain_);
        this->additional_matrix[i] = std::move(row_vector);
    }

    // compute the encoding of additional matrix
    this->matrix_add_commit.resize(this->codeword_domain_size_);
    for(std::size_t i=0; i<this->codeword_domain_size_; i++){
        std::vector<FieldT> add_slice(this->RS_col_systematic_domain_size_);
        for(std::size_t j=0;j<this->RS_col_systematic_domain_size_;j++){
            add_slice[j]=this->additional_matrix[j][i];
        }
        const std::vector<FieldT> add_col_coefficients=IFFT_over_field_subset<FieldT>(add_slice, this->RS_col_systematic_domain_);
        std::vector<FieldT> add_col_vector=FFT_over_field_subset(add_col_coefficients, this->RS_col_codeword_domain_);
        this->matrix_add_commit[i]=std::move(add_col_vector);
    }

    libff::leave_block("Calculating encoding matrix of additional_input");

    libff::enter_block("Calculating and submitting response: Lincheck A");


    this->lincheck_A_->calculate_and_submit_responses(this->additional_matrix,
                                                     this->lincheck_ABC_randomcombinations,
                                                     this->rowcheck_blinding_vector,
                                                     this->query_set,
                                                     this->Uw_matrix,
                                                     this->Ux_matrix);
    libff::leave_block("Calculating and submitting response: Lincheck A");

    libff::enter_block("Calculating and submitting response: Lincheck B");
    this->lincheck_B_->calculate_and_submit_responses(this->additional_matrix,
                                                     this->lincheck_ABC_randomcombinations,
                                                      this->rowcheck_blinding_vector,
                                                     this->query_set,
                                                     this->Uw_matrix,
                                                     this->Uy_matrix);
    libff::leave_block("Calculating and submitting response: Lincheck B");

    libff::enter_block("Calculating and submitting response: Lincheck C");
    this->lincheck_C_->calculate_and_submit_responses(this->additional_matrix,
                                                     this->lincheck_ABC_randomcombinations,
                                                      this->rowcheck_blinding_vector,
                                                     this->query_set,
                                                     this->Uw_matrix,
                                                     this->Uz_matrix);
    libff::leave_block("Calculating and submitting response: Lincheck C");
//编码后矩阵按列存储
    libff::enter_block("Calculating and submitting response: IRScheck");
    this->IRScheck_w_->calculate_and_submit_responses(this->Uw_matrix,
                                                      this->query_set,
                                                      this->IRS_random_combinations,
                                                      this->rowcheck_blinding_matrix);

    this->IRScheck_x_->calculate_and_submit_responses(this->Ux_matrix,
                                                      this->query_set,
                                                      this->IRS_random_combinations,
                                                      this->IRScheck_w_->public_polys,
                                                      this->rowcheck_blinding_matrix);

    this->IRScheck_y_->calculate_and_submit_responses(this->Uy_matrix,
                                                      this->query_set,
                                                      this->IRS_random_combinations,
                                                      this->IRScheck_w_->public_polys,
                                                      this->rowcheck_blinding_matrix);

    this->IRScheck_z_->calculate_and_submit_responses(this->Uz_matrix,
                                                      this->query_set,
                                                      this->IRS_random_combinations,
                                                      this->IRScheck_w_->public_polys,
                                                      this->rowcheck_blinding_matrix);

    libff::leave_block("Calculating and submitting response: IRScheck");

    libff::enter_block("Calculating and submitting response: Rowcheck");

    this->rowcheck_->calculate_and_submit_responses(this->Ux_matrix,
                                                    this->Uy_matrix,
                                                    this->Uz_matrix,
                                                    this->query_set,
                                                    this->IRS_random_combinations,
                                                    this->IRScheck_w_->public_polys,
                                                    this->rowcheck_blinding_vector);

    std::vector<std::vector<FieldT>> ().swap(this->Uw_matrix);
    std::vector<std::vector<FieldT>> ().swap(this->Ux_matrix);
    std::vector<std::vector<FieldT>> ().swap(this->Uy_matrix);
    std::vector<std::vector<FieldT>> ().swap(this->Uz_matrix);

    libff::leave_block("Calculating and submitting response: Rowcheck");
    float ipa_time = inner_product_prover(this->parameters_);
    libff::leave_block("Calculating and submitting proof");
    return (random_time+ipa_time);
}


/* Verification */
template<typename FieldT>
bool interleaved_r1cs_protocol<FieldT>::verifier_predicate(const std::map<std::size_t,FieldT> &public_input)
{
    libff::enter_block("commitment check");

    bool suc=this->merkleTree->verify_merkle_commit(this->lincheck_blinding_commit[0]);
    if(!(suc)){
        std::cout<<"false\n";
        return false;
    }

    libff::leave_block("commitment check");

    std::vector<FieldT> additional_input = std::vector<FieldT>(this->matrix_width_, FieldT(0));//w_bar
    typename std::map<std::size_t, FieldT>::const_iterator it;
    for (it = public_input.begin(); it != public_input.end(); it++)
    {
        std::size_t idx = it->first;
        FieldT val = it->second;
	    additional_input[idx] = val;
	//else public_input has sth wrong?
    }//v_bar
    const std::size_t additional_input_size = additional_input.size();
    const std::size_t target_size = this->num_oracles_vectors_ * this->systematic_domain_size_;
    std::vector<FieldT> additional_target = std::vector<FieldT>(target_size, FieldT(0));

    libff::enter_block("Checking predicate for Lincheck A");
    if (!this->lincheck_A_->verifier_predicate(additional_input,
                                               this->lincheck_ABC_randomcombinations,
                                               this->lincheck_A_->response_polys_coefficients,
                                               this->query_set))
    {
        libff::print_indent(); printf("Interleaved Lincheck for A matrix failed\n");
        return false;
    }
    libff::leave_block("Checking predicate for Lincheck A");

    libff::enter_block("Checking predicate for Lincheck B");
    if (!this->lincheck_B_->verifier_predicate(additional_input,
                                               this->lincheck_ABC_randomcombinations,
                                               this->lincheck_B_->response_polys_coefficients,
                                               this->query_set))
    {
        libff::print_indent(); printf("Interleaved Lincheck for B matrix failed\n");
        return false;
    }
    libff::leave_block("Checking predicate for Lincheck B");

    libff::enter_block("Checking predicate for Lincheck C");
    if (!this->lincheck_C_->verifier_predicate(additional_input,
                                               this->lincheck_ABC_randomcombinations,
                                               this->lincheck_C_->response_polys_coefficients,
                                               this->query_set))
    {
        libff::print_indent(); printf("Interleaved Lincheck for C matrix failed\n");
        return false;
    }
    libff::leave_block("Checking predicate for Lincheck C");

    libff::enter_block("Checking predicate for IRScheck");
    if(!this->IRScheck_w_->verifier_predicate(this->IRS_random_combinations,
                                              this->IRScheck_w_->evals_of_response_polys,
                                              this->query_set)){
        libff::print_indent();
        printf("IRScheck Uw failed\n");
        return false;
    }
    if(!this->IRScheck_x_->verifier_predicate(this->IRS_random_combinations,
                                              this->IRScheck_x_->evals_of_response_polys,
                                              this->query_set)){
        libff::print_indent();
        printf("IRScheck Ux failed\n");
        return false;
    }
    if(!this->IRScheck_y_->verifier_predicate(this->IRS_random_combinations,
                                              this->IRScheck_y_->evals_of_response_polys,
                                              this->query_set)){
        libff::print_indent();
        printf("IRScheck Uy failed\n");
        return false;
    }
    if(!this->IRScheck_z_->verifier_predicate(this->IRS_random_combinations,
                                              this->IRScheck_z_->evals_of_response_polys,
                                              this->query_set)){
        libff::print_indent();
        printf("IRScheck Uz failed\n");
        return false;
    }
    libff::leave_block("Checking predicate for IRScheck");

    libff::enter_block("Checking predicate for Rowcheck");
    if (!this->rowcheck_->verifier_predicate(this->rowcheck_randomcombinations,
                                             this->rowcheck_->response_coefficients,
                                             this->query_set))
    {
        libff::print_indent(); printf("Interleaved Rowcheck failed\n");
        return false;
    }
    libff::leave_block("Checking predicate for Rowcheck");

    if(!inner_product_verifier(this->parameters_)){
        return false;
    }
    return true;
}

template<typename FieldT>
float interleaved_r1cs_protocol<FieldT>::inner_product_prover(const encoded_ligero_parameters &parameters){
    // TODO OPT: 可以有更好的方法实现向量拼接 减少变量的复制

    struct timeval start1,end1;
    gettimeofday(&start1, nullptr);

    libff::enter_block("initial polys");
    std::vector<polynomial<FieldT>> IPA_sec_polys;
    std::vector<polynomial<FieldT>> IPA_pub_polys;
    FieldT target_sum(0);

    IPA_sec_polys.insert(IPA_sec_polys.end(),this->IRScheck_w_->secret_polys.begin(),this->IRScheck_w_->secret_polys.end());
    IPA_pub_polys.insert(IPA_pub_polys.end(),this->IRScheck_w_->public_polys.begin(),this->IRScheck_w_->public_polys.end());
    target_sum+=this->IRScheck_w_->target_sum;

    std::vector<polynomial<FieldT>> ().swap(this->IRScheck_w_->secret_polys);
    std::vector<polynomial<FieldT>> ().swap(this->IRScheck_w_->public_polys);

    IPA_sec_polys.insert(IPA_sec_polys.end(),this->IRScheck_x_->secret_polys.begin(),this->IRScheck_x_->secret_polys.end());
    IPA_pub_polys.insert(IPA_pub_polys.end(),this->IRScheck_x_->public_polys.begin(),this->IRScheck_x_->public_polys.end());
    target_sum+=this->IRScheck_x_->target_sum;

    std::vector<polynomial<FieldT>> ().swap(this->IRScheck_x_->secret_polys);
    std::vector<polynomial<FieldT>> ().swap(this->IRScheck_x_->public_polys);

    IPA_sec_polys.insert(IPA_sec_polys.end(),this->IRScheck_y_->secret_polys.begin(),this->IRScheck_y_->secret_polys.end());
    IPA_pub_polys.insert(IPA_pub_polys.end(),this->IRScheck_y_->public_polys.begin(),this->IRScheck_y_->public_polys.end());
    target_sum+=this->IRScheck_y_->target_sum;

    std::vector<polynomial<FieldT>> ().swap(this->IRScheck_y_->secret_polys);
    std::vector<polynomial<FieldT>> ().swap(this->IRScheck_y_->public_polys);

    IPA_sec_polys.insert(IPA_sec_polys.end(),this->IRScheck_z_->secret_polys.begin(),this->IRScheck_z_->secret_polys.end());
    IPA_pub_polys.insert(IPA_pub_polys.end(),this->IRScheck_z_->public_polys.begin(),this->IRScheck_z_->public_polys.end());
    target_sum+=this->IRScheck_z_->target_sum;

    std::vector<polynomial<FieldT>> ().swap(this->IRScheck_z_->secret_polys);
    std::vector<polynomial<FieldT>> ().swap(this->IRScheck_z_->public_polys);

    IPA_sec_polys.insert(IPA_sec_polys.end(),this->lincheck_A_->secret_polys.begin(),this->lincheck_A_->secret_polys.end());
    IPA_pub_polys.insert(IPA_pub_polys.end(),this->lincheck_A_->public_polys.begin(),this->lincheck_A_->public_polys.end());
    target_sum+=this->lincheck_A_->target_sum;

    std::vector<polynomial<FieldT>> ().swap(this->lincheck_A_->secret_polys);
    std::vector<polynomial<FieldT>> ().swap(this->lincheck_A_->public_polys);

    IPA_sec_polys.insert(IPA_sec_polys.end(),this->lincheck_B_->secret_polys.begin(),this->lincheck_B_->secret_polys.end());
    IPA_pub_polys.insert(IPA_pub_polys.end(),this->lincheck_B_->public_polys.begin(),this->lincheck_B_->public_polys.end());
    target_sum+=this->lincheck_B_->target_sum;

    std::vector<polynomial<FieldT>> ().swap(this->lincheck_B_->secret_polys);
    std::vector<polynomial<FieldT>> ().swap(this->lincheck_B_->public_polys);

    IPA_sec_polys.insert(IPA_sec_polys.end(),this->lincheck_C_->secret_polys.begin(),this->lincheck_C_->secret_polys.end());
    IPA_pub_polys.insert(IPA_pub_polys.end(),this->lincheck_C_->public_polys.begin(),this->lincheck_C_->public_polys.end());
    target_sum+=this->lincheck_C_->target_sum;

    std::vector<polynomial<FieldT>> ().swap(this->lincheck_C_->secret_polys);
    std::vector<polynomial<FieldT>> ().swap(this->lincheck_C_->public_polys);

    IPA_sec_polys.insert(IPA_sec_polys.end(),this->rowcheck_->secret_polys.begin(),this->rowcheck_->secret_polys.end());
    IPA_pub_polys.insert(IPA_pub_polys.end(),this->rowcheck_->public_polys.begin(),this->rowcheck_->public_polys.end());
    target_sum+=this->rowcheck_->target_sum;

    std::vector<polynomial<FieldT>> ().swap(this->rowcheck_->secret_polys);
    std::vector<polynomial<FieldT>> ().swap(this->rowcheck_->public_polys);

    libff::leave_block("initial polys");

    gettimeofday(&end1, nullptr);
    const float initial_time = (end1.tv_usec-start1.tv_usec)/1000000.0 + end1.tv_sec-start1.tv_sec;

    libff::enter_block("Setting parameters");
    // TODO 用下面的代码检查内积关系是否正确 可能多项式的度大于求值域 需要换一个方法
//    FieldT val1(0);
//    for(std::size_t i=0;i<IPA_sec_polys.size();i++){
//        std::vector<FieldT> secret_eval=IPA_sec_polys[i].evaluations_over_field_subset(this->RS_col_systematic_domain_);
//        std::vector<FieldT> public_eval=IPA_pub_polys[i].evaluations_over_field_subset(this->RS_col_systematic_domain_);
//        for(std::size_t j=0;j<secret_eval.size();j++){
//            val1+=secret_eval[j]*public_eval[j];
//        }
//    }
//    assert(val1==target_sum);

    // TODO 这里可能不是相加 而应该求最大值
    std::size_t summation_degree_bound=IPA_pub_polys[0].degree()+IPA_sec_polys[0].degree()+1;
    std::size_t poly_degree_bound=libff::round_to_next_power_of_2(IPA_pub_polys[0].degree()+IPA_sec_polys[0].degree());
    std::vector<std::size_t> localization_parameter_array = this->parameters_.localization_parameter_array;
    //localization_parameter_array.insert(localization_parameter_array.begin(), 1);

    // the summation domain, inner product domain
    const std::size_t compute_domain_dim = libff::log2(this->RS_col_systematic_domain_.num_elements());

    // after encoding
    const std::size_t IPA_codeword_domain_dim = compute_domain_dim + IPA_RS_extra_dimensions_;
    /** NOTE: must be codeword_domain (1<<size, FieldT(1<<size))
     * do not know why now...**/
    // vertical codeword domain
    field_subset<FieldT> IPA_codeword_domain(1 << IPA_codeword_domain_dim, FieldT(1 << IPA_codeword_domain_dim));
    const std::size_t codeword_domain_size = 1 << IPA_codeword_domain_dim;

    const long double field_size_bits = (long double)(libff::soundness_log_of_field_size_helper<FieldT>(FieldT::zero()));
    const std::size_t inter_repetition_num = this->parameters_.IPA_num_interaction_repetitions_;
    const std::size_t query_repetition_num = this->parameters_.IPA_num_query_repetitions_;
    IPA_query_set.clear();
    for (std::size_t i = 0; i < query_repetition_num; ++i)
    {
        bool is_repeat = true;
        std::size_t val;
        while(is_repeat){
            val = std::rand() % (codeword_domain_size);
            //val = std::rand() % (codeword_domain_size >> localization_parameter_array[0]);
            std::vector<std::size_t>::iterator it;
            it = find(IPA_query_set.begin(),IPA_query_set.end(),val);
            if(it == IPA_query_set.end()){
                is_repeat = false;
            }
        }
        IPA_query_set.emplace_back(val);
    }
    libff::leave_block("Setting parameters");

    libff::enter_block("FFTs for pub polys");

    // TODO FFT delegation
    // NOTE: these FFTs should be computed by the verifier, here I have adjusted to let the prover compute this time
    // The only remained to do is let the prover prove to the verifier the FFT values are right
    std::vector<std::vector<FieldT>> IPA_pub_evaluations;

    IPA_pub_evaluations.resize(IPA_pub_polys.size());
    std::size_t set_size = (IPA_pub_polys.size()>>4);

    float p_fft_time = 0;
    float p_fft_GKR_prove_time = 0;
    size_t fft_GKR_proof_size_field = 0;

    std::vector<FieldT> poly_coe_GKR(IPA_pub_polys[0].coefficients().size(),FieldT::zero());
    std::vector<FieldT> poly_eva_GKR(this->RS_col_codeword_domain_size_, FieldT::zero());

    struct timeval start3,end3;
    gettimeofday(&start3, nullptr);

    for (std::size_t i = 0; i < (set_size); i ++)
    {
        IPA_pub_evaluations[i] = FFT_over_field_subset(IPA_pub_polys[i].coefficients(),this->RS_col_codeword_domain_);
        IPA_pub_evaluations[i+ set_size] = FFT_over_field_subset(IPA_pub_polys[i+ set_size].coefficients(),this->RS_col_codeword_domain_);

        IPA_pub_evaluations[i + set_size*2] = IPA_pub_evaluations[i];
        IPA_pub_evaluations[i + set_size*3] = IPA_pub_evaluations[i+ set_size];
        IPA_pub_evaluations[i + set_size*4] = IPA_pub_evaluations[i];
        IPA_pub_evaluations[i + set_size*5] = IPA_pub_evaluations[i+ set_size];
        IPA_pub_evaluations[i + set_size*6] = IPA_pub_evaluations[i];
        IPA_pub_evaluations[i + set_size*7] = IPA_pub_evaluations[i+ set_size];
        IPA_pub_evaluations[i + set_size*14] = IPA_pub_evaluations[i];
        IPA_pub_evaluations[i + set_size*15] = IPA_pub_evaluations[i+ set_size];
    }

    for (std::size_t i = (set_size * 8); i < (set_size * 14); i ++)
    {
        IPA_pub_evaluations[i] = FFT_over_field_subset(IPA_pub_polys[i].coefficients(),this->RS_col_codeword_domain_);
    }

    gettimeofday(&end3, nullptr);
    const float FFT_pub_time = (end3.tv_usec-start3.tv_usec)/1000000.0 + end3.tv_sec-start3.tv_sec ;
    std::cout << "FFT_pub_time is " << FFT_pub_time << std::endl;

    struct timeval start2,end2;
    gettimeofday(&start2, nullptr);

    for (std::size_t j = 0; j < IPA_pub_polys[0].coefficients().size(); j ++)
    {
        for (std::size_t i = 0; i < (set_size)*16; i ++)
        {
            poly_coe_GKR[j] += IPA_pub_polys[i].coefficients()[j];
        }
    }

    for (std::size_t j = 0; j < IPA_pub_evaluations[0].size(); j ++)
    {
        for (std::size_t i = 0; i < (set_size)*16; i ++)
        {
            poly_eva_GKR[j] += IPA_pub_evaluations[i][j];
        }
    }

    fft_GKRs.emplace_back(this->RS_col_codeword_domain_.coset(),0);

    fft_GKRs[0].prover_compute(poly_coe_GKR,poly_eva_GKR);
    p_fft_time += fft_GKRs[0].p_fft_time;
    p_fft_GKR_prove_time += fft_GKRs[0].p_time;
    fft_GKR_proof_size_field += fft_GKRs[0].proof_size;
//
//    proof size compute for 64 and 128
//    p_fft_time += 0;
//    p_fft_GKR_prove_time += 0;
//    fft_GKR_proof_size_field += 10;

    gettimeofday(&end2, nullptr);
    const float setting_GKR_time = (end2.tv_usec-start2.tv_usec)/1000000.0 + end2.tv_sec-start2.tv_sec - p_fft_GKR_prove_time;

    libff::leave_block("FFTs for pub polys");

    // TODO FFT delegation

    struct timeval start,end;
    gettimeofday(&start, nullptr);

    libff::enter_block("FFTs for sec polys");
    // The FFT evluation for secret poly
//    std::vector<std::vector<FieldT>> IPA_sec_polys_evluation_on_codeword_domain;
//    IPA_sec_polys_evluation_on_codeword_domain.resize(IPA_sec_polys.size());
//    for (std::size_t i = 0; i < IPA_sec_polys.size() ; i++) {
//        IPA_sec_polys_evluation_on_codeword_domain[i] = (FFT_over_field_subset(IPA_sec_polys[i].coefficients(), this->RS_col_codeword_domain_));
//    }

    // TODO substitute the secret vectors evaluation and send it to the verifier
    // TODO ldt verifier improve
    // TODO linecheck matrix
    // TODO repetition of public and secret matrices

    std::vector<std::vector<FieldT>> IPA_sec_evaluations;
    IPA_sec_evaluations.resize(IPA_sec_polys.size());
    // 512
    assert( set_size = (IPA_pub_polys.size() >> 4));
    std::size_t size_2 = (IPA_sec_polys.size()>>4) * 2;
    std::size_t size_3 = (IPA_sec_polys.size()>>4) * 3;
    std::size_t size_4 = (IPA_sec_polys.size()>>4) * 4;
    std::size_t size_5 = (IPA_sec_polys.size()>>4) * 5;
    std::size_t size_6 = (IPA_sec_polys.size()>>4) * 6;
    std::size_t size_7 = (IPA_sec_polys.size()>>4) * 7;
    std::size_t size_8 = (IPA_sec_polys.size()>>4) * 8;
    std::size_t size_9 = (IPA_sec_polys.size()>>4) * 9;
    std::size_t size_10 = (IPA_sec_polys.size()>>4) * 10;
    std::size_t size_11 = (IPA_sec_polys.size()>>4) * 11;
    std::size_t size_12 = (IPA_sec_polys.size()>>4) * 12;
    std::size_t size_13 = (IPA_sec_polys.size()>>4) * 13;
    std::size_t size_14 = (IPA_sec_polys.size()>>4) * 14;
    std::size_t size_15 = (IPA_sec_polys.size()>>4) * 15;


    for (std::size_t i = 0; i< set_size; i ++)
    {
        IPA_sec_evaluations[i].resize(this->matrix_w_commit[0].size());
        IPA_sec_evaluations[i+size_2].resize(this->matrix_w_commit[0].size());
        IPA_sec_evaluations[i+size_4].resize(this->matrix_w_commit[0].size());
        IPA_sec_evaluations[i+size_6].resize(this->matrix_w_commit[0].size());
        IPA_sec_evaluations[i+size_8].resize(this->matrix_w_commit[0].size());
        IPA_sec_evaluations[i+size_10].resize(this->matrix_w_commit[0].size());
        IPA_sec_evaluations[i+size_12].resize(this->matrix_w_commit[0].size());

        for (std::size_t j = 0; j < this->matrix_w_commit[0].size(); j ++)
        {
            IPA_sec_evaluations[i][j] = (this->matrix_w_commit[i][j] + matrix_Uxyz_Rowcheck_blinding_commit[i][j]);
            //test_evaluations[i+ set_size][j] = test_evaluations[i][j];
            IPA_sec_evaluations[i+ size_2][j] = (this->matrix_x_commit[i][j] + matrix_Uxyz_Rowcheck_blinding_commit[i][j]);
            //test_evaluations[i+ size_3][j] = test_evaluations[i+ size_2][j];
            IPA_sec_evaluations[i+ size_4][j] = (this->matrix_y_commit[i][j] + matrix_Uxyz_Rowcheck_blinding_commit[i][j]);
            //test_evaluations[i+ size_5][j] = test_evaluations[i+ size_4][j];
            IPA_sec_evaluations[i+ size_6][j] = (this->matrix_z_commit[i][j] + matrix_Uxyz_Rowcheck_blinding_commit[i][j]);
//            test_evaluations[i+ size_7][j] = test_evaluations[i+ size_6][j];
            IPA_sec_evaluations[i+ size_8][j] = (this->matrix_w_commit[i][j] + this->matrix_x_commit[i][j] + matrix_Uxyz_Rowcheck_blinding_commit[i][j]
                    + this->matrix_add_commit[i][j]);
            //test_evaluations[i+ size_9][j] = test_evaluations[i+ size_8][j];
            IPA_sec_evaluations[i+ size_10][j] = (this->matrix_w_commit[i][j] + this->matrix_y_commit[i][j] + matrix_Uxyz_Rowcheck_blinding_commit[i][j]
                                                  + this->matrix_add_commit[i][j]);
            //test_evaluations[i+ size_11][j] = test_evaluations[i+ size_10][j];
            IPA_sec_evaluations[i+ size_12][j] = (this->matrix_w_commit[i][j] + this->matrix_z_commit[i][j] + matrix_Uxyz_Rowcheck_blinding_commit[i][j]
                                                   + this->matrix_add_commit[i][j]);
//            test_evaluations[i+ size_13][j] = test_evaluations[i+ size_12][j];
        }
        IPA_sec_evaluations[i+ set_size] = IPA_sec_evaluations[i];
        IPA_sec_evaluations[i+ size_3] = IPA_sec_evaluations[i+ size_2];
        IPA_sec_evaluations[i+ size_5] = IPA_sec_evaluations[i+ size_4];
        IPA_sec_evaluations[i+ size_7] = IPA_sec_evaluations[i+ size_6];
        IPA_sec_evaluations[i+ size_9] = IPA_sec_evaluations[i+ size_8];
        IPA_sec_evaluations[i+ size_11] = IPA_sec_evaluations[i+ size_10];
        IPA_sec_evaluations[i+ size_13] = IPA_sec_evaluations[i+ size_12];
        IPA_sec_evaluations[i+ size_14] = (FFT_over_field_subset(IPA_sec_polys[i+ size_14].coefficients(), this->RS_col_codeword_domain_));
        IPA_sec_evaluations[i+ size_15] = IPA_sec_evaluations[i+ size_14] ;
    }

    std::vector<std::vector<FieldT>> ().swap(this->matrix_w_commit);
    std::vector<std::vector<FieldT>> ().swap(this->matrix_x_commit);
    std::vector<std::vector<FieldT>> ().swap(this->matrix_y_commit);
    std::vector<std::vector<FieldT>> ().swap(this->matrix_z_commit);
    std::vector<std::vector<FieldT>> ().swap(this->matrix_Uxyz_Rowcheck_blinding_commit);

    libff::leave_block("FFTs for sec polys");

    std::size_t padding_degree = poly_degree_bound - summation_degree_bound;

    libff::enter_block("Setting Inner Product Verifier");

//    this->IPA_verifier_.reset(new Inner_product_verifier<FieldT>(IPA_pub_polys, this->RS_col_systematic_domain_, IPA_sec_evaluations,
//                                                                 IPA_pub_evaluations, padding_degree, poly_degree_bound,
//                                                                 localization_parameter_array, IPA_codeword_domain, target_sum, inter_repetition_num));

    this->IPA_verifier_.reset(new Inner_product_verifier<FieldT>(std::move(IPA_pub_polys), this->RS_col_systematic_domain_, std::move(IPA_sec_evaluations),
                                                                 std::move(IPA_pub_evaluations), padding_degree, poly_degree_bound,
                                                                 localization_parameter_array, IPA_codeword_domain, target_sum, inter_repetition_num));
    std::vector<std::vector<FieldT>>().swap(IPA_sec_evaluations);
    std::vector<std::vector<FieldT>>().swap(IPA_pub_evaluations);

    libff::leave_block("Setting Inner Product Verifier");

    gettimeofday(&end, nullptr);
    const float setting_time = (end.tv_usec-start.tv_usec)/1000000.0 + end.tv_sec-start.tv_sec + initial_time + setting_GKR_time;

    libff::enter_block("Inner Product Prover");
    libff::enter_block("Setting Inner Product Prover and compute the first round");
    // TODO: There is no need to commit v_trees in the prover
    // TODO: The commitment for evaluation can also be split

    this->IPA_prover_.reset(new Inner_product_prover<FieldT>(std::move(IPA_pub_polys), std::move(IPA_sec_polys), IPA_query_set,
                                                             localization_parameter_array, poly_degree_bound, *this->IPA_verifier_, IPA_codeword_domain,
                                                             inter_repetition_num));

    libff::leave_block("Setting Inner Product Prover and compute the first round");
    libff::enter_block("Proving all the remained rounds for FRI");
    this->IPA_prover_->prove(IPA_query_set);
    libff::leave_block("Proving all the remained rounds for FRI");

    libff::enter_block("Compute proof size");
    // 一共IPA_sec_polys.size个矩阵 每一列的路径长为merklenth
    std::size_t merklenth = this->lincheck_blinding_commit[0].path_lenth ;
    std::size_t v_trees_hashes = 1 * merklenth;
    std::size_t h_tree_hashes = this->IPA_prover_->h_tree_lenth;
    std::size_t FRI_trees_hashes = this->IPA_prover_->FRI_tree_lenth ;
    // could construct H = (h_1, h_2 ,..., h_{codeword})
    std::size_t proof_size_roots_hash_number = 1 + 1 + (localization_parameter_array.size() - 2);
    std::size_t proof_size_path_hash_number = v_trees_hashes  + h_tree_hashes + FRI_trees_hashes;
    std::size_t proof_size_hash_number = proof_size_roots_hash_number + proof_size_path_hash_number;

    std::size_t hash_size = 200;
    std::size_t proof_size_hash = (proof_size_hash_number * hash_size)/1024/8;

    std::size_t proof_size_field_number = 0;
    std::size_t poly_degree_bound_last_round = poly_degree_bound ;
    // Note that it is necessary to query (query number) times in all rounds instead of each round
    for (std::size_t i = 0 ; i < localization_parameter_array.size() ; i++)
    {
        if (i==0)
        {
            // v field number
            proof_size_field_number +=  5 * this->codeword_domain_size_ * this->parameters_.IPA_num_query_repetitions_
                    * ( 1 << localization_parameter_array[i]);
            // delta + h field number
            proof_size_field_number += 2 *  this->parameters_.IPA_num_query_repetitions_ * ( 1 << localization_parameter_array[i]);
            poly_degree_bound_last_round /= ( 1ull << localization_parameter_array[i]);
        }
        else{
            proof_size_field_number +=  this->parameters_.IPA_num_query_repetitions_
                    * ( 1 << localization_parameter_array[i]);
            poly_degree_bound_last_round /= ( 1ull << localization_parameter_array[i]);
        }
    }

    // The final_poly proof size
    std::size_t FRI_field_number_last_round = inter_repetition_num * (poly_degree_bound_last_round + 1) ;
    proof_size_field_number +=  FRI_field_number_last_round;
    std::size_t ipa_field_number = proof_size_field_number;

    // correctness check field proof size
    std::size_t correct_check_field_number = this->num_interactions_ *
            (4 * this->systematic_domain_size_ + 4 * (2*this->systematic_domain_size_-1)) ;
    proof_size_field_number += correct_check_field_number;

    float correct_check_field_size = correct_check_field_number * field_size_bits / 1024 /8;
    float ipa_field_size = ipa_field_number * field_size_bits / 1024 /8;
    //float proof_size_field = (proof_size_field_number * field_size_bits)/1024/8;

    // proof size of GKR protocol itself
    float fft_GKR_proof_size = fft_GKR_proof_size_field/ 1024 /8; //mark
    std::cout << "fft_GKR_proof_size is " << fft_GKR_proof_size << std::endl;
    // proof size of
//    float fft_GKR_output = (this->parameters_.IPA_num_query_repetitions_ * 8 * 16 *
//            (1ull << this->parameters_.localization_parameter_array[0]) * field_size_bits) / 1024 /8; //mark
    // Here 5 means 1 IRS + 3 Linecheck + 1 Rowcheck
    float fft_GKR_output = (this->parameters_.IPA_num_query_repetitions_ * this->parameters_.IPA_num_interaction_repetitions_ * 5 *
                            (1ull << this->parameters_.localization_parameter_array[0]) * field_size_bits) / 1024 /8; //mark
    std::cout << "fft_GKR_output is " << fft_GKR_output << std::endl;



    float proof_size_field = (proof_size_field_number * field_size_bits)/1024/8 ; //mark


    const float proof_size = proof_size_field + proof_size_hash ;
    std::cout << "IPA_sec_polys.size() is " << IPA_sec_polys.size() << std::endl;
//    std::cout << "inner product: proof_size_field is " << ipa_field_size << std::endl;
//    std::cout << "inner product: last poly size is " << FRI_field_number_last_round* field_size_bits/1024/8 << std::endl;
//    std::cout << "correct_check: proof_size_field is " << correct_check_field_size << std::endl;
    std::cout << "ligero++: total proof size field is " << proof_size_field << std::endl;
//    std::cout << "inner product: v_trees hashes is " << v_trees_hashes * 200 / 1024 / 8 << std::endl;
//    std::cout << "inner product: h_tree hash is " << h_tree_hashes * 200 / 1024 / 8 << std::endl;
//    std::cout << "inner product: FRI_trees hash is " << FRI_trees_hashes * 200 / 1024 / 8 << std::endl;
//    std::cout << "inner product: proof_size_roots_hashes is " << proof_size_roots_hash_number * 200 / 1024 / 8 << std::endl;
//    std::cout << "inner product: proof_size_path_hashes is " << proof_size_path_hash_number * 200 / 1024 / 8 << std::endl;
    std::cout << "inner product: proof_size_hash is " << proof_size_hash << std::endl;
    std::cout << "inner product: proof_size_fft_GKR is " << fft_GKR_proof_size + fft_GKR_output << std::endl;// mark
    std::cout<< "ligero++: total_proof_size is "<<proof_size+ fft_GKR_proof_size + fft_GKR_output<<std::endl;
    libff::leave_block("Compute proof size");

    std::cout << "inner product: the time of calculate ffts:" << p_fft_time << std::endl; // mark
    std::cout << "inner product: the time of calculate fft_GKR proof:" << p_fft_GKR_prove_time << std::endl;// mark
    libff::leave_block("Inner Product Prover");
    return setting_time;
}

template<typename FieldT>
bool interleaved_r1cs_protocol<FieldT>::inner_product_verifier(const encoded_ligero_parameters &parameters) {
    libff::enter_block("Inner Product Verifier");
    bool result = this->IPA_verifier_->verify(IPA_query_set, &(*this->IPA_prover_));
    //mark
    float verify_fft_GKR_time = 0;
    for(size_t i; i < fft_GKRs.size(); i++)
    {
        result &= fft_GKRs[i].verifier_predicate();
        verify_fft_GKR_time += fft_GKRs[i].v_time;
    }
    std::cout << "Inner Product Verifier: the time of verify fft_GKR:" <<verify_fft_GKR_time << std::endl;
    if (!result){
        libff::print_indent(); printf("error occurs! \n");
    }
    libff::leave_block("Inner Product Verifier");
    return result;
}

} // namespace ligero

