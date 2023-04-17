#include <cmath>
#include <stdexcept>
#include <libff/common/profiling.hpp>
#include "ligero/algebra/fft.hpp"
#include "ligero/algebra/lagrange.hpp"

namespace ligero{
template<typename FieldT>
interleaved_RScode<FieldT>::interleaved_RScode(
        iop_protocol<FieldT> &IOP,
        const domain_handle &codeword_domain_handle,
        const domain_handle &systematic_domain_handle,
        const domain_handle &RS_col_systematic_domain_handle,
        const std::size_t num_oracles,
        const std::size_t num_queries,//打开的列
        const std::size_t num_interactions,//交互次数
        const bool make_zk,
        const field_subset_type domain_type) :
        IOP_(IOP),
        codeword_domain_handle_(codeword_domain_handle),
        systematic_domain_handle_(systematic_domain_handle),
        RS_col_systematic_domain_handle_(RS_col_systematic_domain_handle),
        num_oracles_(num_oracles),
        num_queries_(num_queries),
        num_interactions_(num_interactions),
        make_zk_(make_zk),
        field_subset_type_(domain_type)
{
    this->codeword_domain_ = this->IOP_.get_domain(this->codeword_domain_handle_);
    this->codeword_domain_size_ = this->codeword_domain_.num_elements();
    this->systematic_domain_ = this->IOP_.get_domain(this->systematic_domain_handle_);
    this->systematic_domain_size_ = this->systematic_domain_.num_elements();
    this->RS_col_systematic_domain_=this->IOP_.get_domain(this->RS_col_systematic_domain_handle_);
    this->RS_col_systematic_domain_size_=this->RS_col_systematic_domain_.num_elements();
    this->response_size_ = this->codeword_domain_size_;
}

//proving
// 为啥这里给的是blinding矩阵-->因为Uw与Uxyz的blinding 所在位置不同 -->OPT: 可以给blinding的vector而不给矩阵,但需要标识blinding应该填到哪一行
template<typename FieldT>
void interleaved_RScode<FieldT>::calculate_and_submit_responses(const std::vector<std::vector<FieldT>> &IRS_matrix,
                                                                const std::vector<std::size_t>& query_set,
                                                                const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                                                const std::vector<std::vector<FieldT>> &blinding_matrix)
{
    //IRS_matrix是按行存储的
    // 随机挑战的长度为RS_sys_domain_size
    this->evals_of_response_polys.resize(this->num_interactions_);
    this->query_col.resize(this->num_interactions_);
    target_sum=FieldT(0);
    std::size_t count=0;
    secret_polys.resize(this->num_interactions_*this->num_queries_);
    public_polys.resize(this->num_interactions_*this->num_queries_);
    for (std::size_t i = 0; i < this->num_interactions_; ++i)
    {
        const std::vector<FieldT> random_linear_combination = random_linear_combinations[i];
        this->evals_of_response_polys[i].resize(this->codeword_domain_size_,FieldT(0));

        // 计算rU
        if(!this->make_zk_){
            for(std::size_t column_idx=0;column_idx<this->codeword_domain_size_;column_idx++){
                for(std::size_t row_idx=0;row_idx<random_linear_combination.size();row_idx++){
                    this->evals_of_response_polys[i][column_idx]+=(random_linear_combination[row_idx]*(IRS_matrix[row_idx][column_idx]));
                }
            }
            //构造内积
            for(unsigned long pos : query_set){
                std::vector<FieldT> secret_vec;
                secret_vec.resize(this->RS_col_systematic_domain_size_);
                for(std::size_t h=0;h<this->RS_col_systematic_domain_size_;h++){
                    secret_vec[h]=IRS_matrix[h][pos];
                }
                polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                polynomial<FieldT> public_poly(IFFT_over_field_subset(random_linear_combination,this->RS_col_systematic_domain_));
                secret_polys[count]=std::move(secret_poly);
                public_polys[count]=std::move(public_poly);
                target_sum+=this->evals_of_response_polys[i][pos];
                //TODO 不能注释
                count++;
                // 内积检查
//                std::vector<FieldT> secret_eval=secret_polys[count-1].evaluations_over_field_subset(this->RS_col_systematic_domain_);
//                std::vector<FieldT> public_eval=public_polys[count-1].evaluations_over_field_subset(this->RS_col_systematic_domain_);
//                FieldT val(0);
//                for(std::size_t h=0;h<this->RS_col_systematic_domain_size_;h++){
//                    val+=secret_eval[h]*public_eval[h];
//                }
//                assert(val==this->evals_of_response_polys[i][pos]);
            }

        }else{
            // blinding矩阵只有一行有值
            for(std::size_t column_idx=0;column_idx<this->codeword_domain_size_;column_idx++){
                for(std::size_t row_idx=0;row_idx<random_linear_combination.size();row_idx++){
                    this->evals_of_response_polys[i][column_idx]+=(random_linear_combination[row_idx]*(IRS_matrix[row_idx][column_idx]+blinding_matrix[row_idx][column_idx]));
                }
            }
            //构造内积
            for(unsigned long pos : query_set){
                /** In each round we use the same secret_polys **/
                // TODO change
                // TODO each round public_poly also has same
                if (i==0)
                {
                    std::vector<FieldT> secret_vec;
                    secret_vec.resize(this->RS_col_systematic_domain_size_);
                    for(std::size_t h=0;h<this->RS_col_systematic_domain_size_;h++){
                        secret_vec[h]=IRS_matrix[h][pos]+blinding_matrix[h][pos];
                    }
                    polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                    polynomial<FieldT> public_poly(IFFT_over_field_subset(random_linear_combination,this->RS_col_systematic_domain_));
                    secret_polys[count]=std::move(secret_poly);
                    public_polys[count]=std::move(public_poly);
                    target_sum+=this->evals_of_response_polys[i][pos];
                    count++;
                }
                else{
                    secret_polys[count]=secret_polys[count - query_set.size()];
                    polynomial<FieldT> public_poly(IFFT_over_field_subset(random_linear_combination,this->RS_col_systematic_domain_));
                    public_polys[count]=std::move(public_poly);
                    target_sum+=this->evals_of_response_polys[i][pos];
                    count++;
                }
            }
        }
    }
}


template<typename FieldT>
void interleaved_RScode<FieldT>::calculate_and_submit_responses(const std::vector<std::vector<FieldT>> &IRS_matrix,
                                                                const std::vector<std::size_t>& query_set,
                                                                const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                                                const std::vector<polynomial<FieldT>> input_public_polys,
                                                                const std::vector<std::vector<FieldT>> &blinding_matrix)
{
    //IRS_matrix是按行存储的
    // 随机挑战的长度为RS_sys_domain_size
    this->evals_of_response_polys.resize(this->num_interactions_);
    this->query_col.resize(this->num_interactions_);
    target_sum=FieldT(0);
    std::size_t count=0;
    secret_polys.resize(this->num_interactions_*this->num_queries_);
    public_polys.resize(this->num_interactions_*this->num_queries_);
    for (std::size_t i = 0; i < this->num_interactions_; ++i)
    {
        const std::vector<FieldT> random_linear_combination = random_linear_combinations[i];
        this->evals_of_response_polys[i].resize(this->codeword_domain_size_,FieldT(0));

        // 计算rU
        if(!this->make_zk_){
            for(std::size_t column_idx=0;column_idx<this->codeword_domain_size_;column_idx++){
                for(std::size_t row_idx=0;row_idx<random_linear_combination.size();row_idx++){
                    this->evals_of_response_polys[i][column_idx]+=(random_linear_combination[row_idx]*(IRS_matrix[row_idx][column_idx]));
                }
            }
            //构造内积
            for(unsigned long pos : query_set){
                std::vector<FieldT> secret_vec;
                secret_vec.resize(this->RS_col_systematic_domain_size_);
                for(std::size_t h=0;h<this->RS_col_systematic_domain_size_;h++){
                    secret_vec[h]=IRS_matrix[h][pos];
                }
                polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                secret_polys[count]=std::move(secret_poly);
                public_polys[count]=std::move(input_public_polys[count]);
                target_sum+=this->evals_of_response_polys[i][pos];
                //TODO 不能注释
                count++;
                // 内积检查
//                std::vector<FieldT> secret_eval=secret_polys[count-1].evaluations_over_field_subset(this->RS_col_systematic_domain_);
//                std::vector<FieldT> public_eval=public_polys[count-1].evaluations_over_field_subset(this->RS_col_systematic_domain_);
//                FieldT val(0);
//                for(std::size_t h=0;h<this->RS_col_systematic_domain_size_;h++){
//                    val+=secret_eval[h]*public_eval[h];
//                }
//                assert(val==this->evals_of_response_polys[i][pos]);
            }

        }else{
            // blinding矩阵只有一行有值
            for(std::size_t column_idx=0;column_idx<this->codeword_domain_size_;column_idx++){
                for(std::size_t row_idx=0;row_idx<random_linear_combination.size();row_idx++){
                    this->evals_of_response_polys[i][column_idx]+=(random_linear_combination[row_idx]*(IRS_matrix[row_idx][column_idx]+blinding_matrix[row_idx][column_idx]));
                }
            }
            //构造内积
            for(unsigned long pos : query_set){
                /** In each round we use the same secret_polys **/
                // TODO change
                // TODO each round public_poly also has same
                if (i==0)
                {
                    std::vector<FieldT> secret_vec;
                    secret_vec.resize(this->RS_col_systematic_domain_size_);
                    for(std::size_t h=0;h<this->RS_col_systematic_domain_size_;h++){
                        secret_vec[h]=IRS_matrix[h][pos]+blinding_matrix[h][pos];
                    }
                    polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                    secret_polys[count]=std::move(secret_poly);
                    public_polys[count]=std::move(input_public_polys[count]);
                    target_sum+=this->evals_of_response_polys[i][pos];
                    count++;
                }
                else{
                    secret_polys[count]=secret_polys[count - query_set.size()];
                    public_polys[count]=std::move(input_public_polys[count]);
                    target_sum+=this->evals_of_response_polys[i][pos];
                    count++;
                }
            }
        }
    }
}

/* Verification */
template<typename FieldT>
bool interleaved_RScode<FieldT>::verifier_predicate(const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                                    const std::vector<std::vector<FieldT>> &response_mes,
                                                    const std::vector<std::size_t> &query_pos)
{
    for (size_t h = 0; h < this->num_interactions_; ++h)
    {
        /* EQUALITY TEST: does the polynomial that was sent sum to the value that it should if the
           claimed statement is true?
           在et中是正确性检查 即判断rAx=rb
           */
        libff::enter_block("IRScheck:correctness check");
        // response vector
        std::vector<FieldT> response = response_mes[h];
        // 试图解码
        const std::vector<FieldT> try_to_decode =
                IFFT_of_known_degree_over_field_subset<FieldT>(response, this->systematic_domain_size_ ,this->codeword_domain_);
        libff::leave_block("IRScheck:correctness check");
    }
    return true;
}
}
