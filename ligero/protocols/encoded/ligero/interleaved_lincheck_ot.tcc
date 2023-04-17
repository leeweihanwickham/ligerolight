#include <cassert>
#include <cmath>
#include <stdexcept>
#include <thread>
#include <libff/common/profiling.hpp>
#include "ligero/algebra/fft.hpp"
#include "interleaved_lincheck_ot.hpp"
#include "ligero/protocols/encoded/ligero/parallel_ot.cpp"
namespace ligero {

template<typename FieldT>
interleaved_lincheck_ot_protocol<FieldT>::interleaved_lincheck_ot_protocol(
        iop_protocol<FieldT> &IOP,
        const domain_handle &codeword_domain_handle,
        const domain_handle &systematic_domain_handle,
        const domain_handle &extended_systematic_domain_handle,
        const domain_handle &RS_col_systematic_domain_handle,
        const std::size_t num_oracles_input,
        const std::size_t num_oracles_target,
        const std::size_t num_queries,
        const std::size_t num_interactions,
        const bool make_zk,
        const field_subset_type domain_type,
        const naive_sparse_matrix<FieldT> constraint_matrix,
        const std::size_t parallel_opt) :
        IOP_(IOP),
        codeword_domain_handle_(codeword_domain_handle),
        systematic_domain_handle_(systematic_domain_handle),
        RS_col_systematic_domain_handle_(RS_col_systematic_domain_handle),
        num_oracles_input_(num_oracles_input),
        num_oracles_target_(num_oracles_target),
        num_queries_(num_queries),
        num_interactions_(num_interactions),
        make_zk_(make_zk),
        field_subset_type_(domain_type),
        constraint_matrix_(constraint_matrix),
        parallel_opt_(parallel_opt)
{
    this->codeword_domain_ = this->IOP_.get_domain(this->codeword_domain_handle_);
    this->codeword_domain_size_ = this->codeword_domain_.num_elements();

    this->systematic_domain_ = this->IOP_.get_domain(this->systematic_domain_handle_);
    this->systematic_domain_size_ = this->systematic_domain_.num_elements();

    this->extended_systematic_domain_ = this->IOP_.get_domain(extended_systematic_domain_handle);

    this->RS_col_systematic_domain_=this->IOP_.get_domain(RS_col_systematic_domain_handle);
    this->RS_col_systematic_domain_size_=this->RS_col_systematic_domain_.num_elements();

    this->response_size_ = 2 * this->systematic_domain_size_;
}


// 计算 r[Pc,-I]*[w+z,x]=0  首先 将Uw与额外输入Uz相加得到U(w+z)，然后再加上Ux
template<typename FieldT>
void interleaved_lincheck_ot_protocol<FieldT>::calculate_and_submit_responses(const std::vector<std::vector<FieldT>>&addition_matrix,//Uz
                                                                              const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                                                              const std::vector<std::vector<FieldT>> &blinding_vectors,
                                                                              const std::vector<std::size_t>& query_set,
                                                                              const std::vector<std::vector<FieldT>> &input_matrix,//Uw
                                                                              const std::vector<std::vector<FieldT>>&target_matrix//Ux
                                                                              )
{


    this->response_polys_coefficients.resize(this->num_interactions_);

    target_sum=FieldT(0);
    std::size_t count=0;
    secret_polys.resize(this->num_interactions_*this->num_queries_);
    public_polys.resize(this->num_interactions_*this->num_queries_);
    for (size_t h = 0; h < this->num_interactions_; ++h)
    {
        const std::vector<FieldT> random_linear_combination = random_linear_combinations[h];
        // \sum s_i*f_i
        std::vector<FieldT> evals_of_response_poly(this->codeword_domain_size_, FieldT(0));

        std::vector<FieldT> row_vector(this->RS_col_systematic_domain_size_ * (this->systematic_domain_size_ ), FieldT(0));// s_1||s_2||...||s_(m_2);

        if(this->make_zk_){
            // Px高度m1*(h-t) 宽度m2(h-t) 填充-I m1*(h-t)
            // 考虑到blinding 首先在Px右边填充宽度为(h-t)的0 然后再填-I
            const std::size_t  index=(this->num_oracles_input_+1) * (this->systematic_domain_size_ );
            // TODO 两种方法 1：不填充Px 直接在s后面添加 <-r> row_vector实际是[s,-r]   2:给Px填充-I
            for (size_t j = 0; j < this->constraint_matrix_.size(); ++j)
            {
                std::map<std::size_t, FieldT> row = this->constraint_matrix_[j];
                typename std::map<std::size_t, FieldT>::iterator it;
                FieldT a;
                for (it = row.begin(); it != row.end(); it++)
                {
                    const std::size_t idx = it->first;
                    const FieldT val = it->second;
                    row_vector[idx] += random_linear_combination[j] * val;
                }
                row_vector[j+index]=-random_linear_combination[j];
            }

            // 为啥不先将Uw Ux blinding  addition_matrix相加? 这样需要设置一个变量 而不是常量const  涉及矩阵的复制  到底是下面的快还是直接复制矩阵快?
            // Uw
            std::vector<std::vector<FieldT>> query_s_matrix;
            query_s_matrix.resize(this->RS_col_systematic_domain_size_,std::vector<FieldT>(this->codeword_domain_size_,FieldT(0)));
            for(std::size_t i=0;i < this->num_oracles_input_;++i){
                const std::size_t start = (this->systematic_domain_size_ ) * i;
                const std::size_t end = (this->systematic_domain_size_ ) + start;

                // row_vector s or r
                const typename std::vector<FieldT>::const_iterator row_first = row_vector.begin() + start;
                const typename std::vector<FieldT>::const_iterator row_last = row_vector.begin() + end;
                std::vector<FieldT> current_vector(row_first, row_last);

                const std::vector<FieldT> poly_coefficients =
                        IFFT_over_field_subset<FieldT>(current_vector, this->systematic_domain_);
                const std::vector<FieldT> current_evaluations =
                        FFT_over_field_subset<FieldT>(poly_coefficients, this->codeword_domain_);
                for (size_t a = 0; a < this->codeword_domain_size_; ++a)
                {
                    evals_of_response_poly[a] += current_evaluations[a] *
                                                 (input_matrix[i][a] + addition_matrix[i][a]);// s()*f_w()
                }
                for(std::size_t k=0;k<this->num_queries_;k++){
                    query_s_matrix[i][k]=current_evaluations[query_set[k]];
                }
            }
            // Uxyz
            for(std::size_t i=this->num_oracles_input_+1;i<this->num_oracles_input_+this->num_oracles_target_+1;++i){
                const std::size_t start = (this->systematic_domain_size_  ) * i;
                const std::size_t end = (this->systematic_domain_size_  ) + start;

                // row_vector s or r
                const typename std::vector<FieldT>::const_iterator row_first = row_vector.begin() + start;
                const typename std::vector<FieldT>::const_iterator row_last = row_vector.begin() + end;
                std::vector<FieldT> current_vector(row_first, row_last);

                const std::vector<FieldT> poly_coefficients =
                        IFFT_over_field_subset<FieldT>(current_vector, this->systematic_domain_);
                const std::vector<FieldT> current_evaluations =
                        FFT_over_field_subset<FieldT>(poly_coefficients, this->codeword_domain_);
                for (size_t a = 0; a < this->codeword_domain_size_; ++a)
                {
                    evals_of_response_poly[a] += current_evaluations[a] * target_matrix[i][a];// s()*f_w()
                }
                for(std::size_t k=0;k<this->num_queries_;k++){
                    query_s_matrix[i][k]=current_evaluations[query_set[k]];
                }
            }
            // blinding
            for(std::size_t i=this->num_oracles_input_;i<this->num_oracles_input_+this->num_oracles_target_+2;i=i+this->num_oracles_target_+1){
                const std::size_t start = (this->systematic_domain_size_ ) * i;
                const std::size_t end = (this->systematic_domain_size_  ) + start;
                const typename std::vector<FieldT>::const_iterator row_first = row_vector.begin() + start;
                const typename std::vector<FieldT>::const_iterator row_last = row_vector.begin() + end;
                std::vector<FieldT> current_vector(row_first, row_last);

                const std::vector<FieldT> poly_coefficients =
                        IFFT_over_field_subset<FieldT>(current_vector, this->systematic_domain_);
                const std::vector<FieldT> current_evaluations =
                        FFT_over_field_subset<FieldT>(poly_coefficients, this->codeword_domain_);
                for (size_t a = 0; a < this->codeword_domain_size_; ++a)
                {
                    evals_of_response_poly[a] += current_evaluations[a] * blinding_vectors[h][a];
                }
                for(std::size_t k=0;k<this->num_queries_;k++){
                    query_s_matrix[i][k]=current_evaluations[query_set[k]];
                }
            }
            // 构造内积
            for(std::size_t k=0;k<this->num_queries_;k++){
                std::size_t pos=query_set[k];
                std::vector<FieldT> public_vec;
                public_vec.resize(this->RS_col_systematic_domain_size_);

                if (h==0){
                    std::vector<FieldT> secret_vec;
                    secret_vec.resize(this->RS_col_systematic_domain_size_);
                    for(std::size_t i=0;i<this->RS_col_systematic_domain_size_;i++){
                        secret_vec[i]=input_matrix[i][pos]+addition_matrix[i][pos]+target_matrix[i][pos];
                        public_vec[i]=query_s_matrix[i][k];
                    }
                    //secret_vec[this->num_oracles_input_]=blinding_vectors[h][pos];
                    secret_vec[this->num_oracles_input_+this->num_oracles_target_+1]=blinding_vectors[h][pos];
                    polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                    polynomial<FieldT> public_poly(IFFT_over_field_subset(public_vec,this->RS_col_systematic_domain_));
                    secret_polys[count]=std::move(secret_poly);
                    public_polys[count]=std::move(public_poly);
                    target_sum+=evals_of_response_poly[pos];
                    count++;
                }
                else{
                    for(std::size_t i=0;i<this->RS_col_systematic_domain_size_;i++){
                        public_vec[i]=query_s_matrix[i][k];
                    }
                    secret_polys[count] = secret_polys[count - this->num_queries_];
                    polynomial<FieldT> public_poly(IFFT_over_field_subset(public_vec,this->RS_col_systematic_domain_));
                    public_polys[count]=std::move(public_poly);
                    target_sum+=evals_of_response_poly[pos];
                    count++;
                }
            }

        } else{
            // TODO No make_zk
            // 算r[Px,-I]=s
            // 首先填充Px 此时的Px: m1h*m2h 填充后宽度应该和Uw的高度一致 因为Uw+Ux的高度和Uw一样
            // Px按行存储，每行是映射，即<index,value>，所以每一行取0的地方不在映射里
            // 根据row_vector的大小推测Px的行长为num_oracles_input_ * systematic_domain_size_ 也就是m2h
            // 在Px右边增加m1h*m1h的负单位阵(行长加m1h)
            // 最后填充Px的行长到RS_col_systematic_domain_size 因为Uw也是填充之后的高度 （这一步可以省略，因为Px不存放值为0的位置）

            // TODO 两种方法 1：不填充Px 直接在s后面添加 <-r>   2:给Px填充-I
            std::size_t  index=this->num_oracles_input_ * (this->systematic_domain_size_ );
            for (size_t j = 0; j < this->constraint_matrix_.size(); ++j)
            {
                std::map<std::size_t, FieldT> row = this->constraint_matrix_[j];
                typename std::map<std::size_t, FieldT>::iterator it;
                FieldT a;
                for (it = row.begin(); it != row.end(); it++)
                {
                    const std::size_t idx = it->first;
                    const FieldT val = it->second;
                    row_vector[idx] += random_linear_combination[j] * val;
                }
                row_vector[j+index]=-random_linear_combination[j];
            }
            // TODO OPT：这里当i<m2时，target_matrix的g每一行是0 可以不算它 时间会减少多少？
            std::vector<std::vector<FieldT>> query_s_matrix;
            query_s_matrix.resize(this->RS_col_systematic_domain_size_,std::vector<FieldT>(this->codeword_domain_size_,FieldT(0)));
            for(std::size_t i=0;i<this->RS_col_systematic_domain_size_;i++){
                const std::size_t start = (this->systematic_domain_size_ ) * i;
                const std::size_t end = (this->systematic_domain_size_) + start;
                // row_vector s or r
                const typename std::vector<FieldT>::const_iterator row_first = row_vector.begin() + start;
                const typename std::vector<FieldT>::const_iterator row_last = row_vector.begin() + end;
                std::vector<FieldT> current_vector(row_first, row_last);
                const std::vector<FieldT> poly_coefficients =
                        IFFT_over_field_subset<FieldT>(current_vector, this->systematic_domain_);
                const std::vector<FieldT> current_evaluations =
                        FFT_over_field_subset<FieldT>(poly_coefficients, this->codeword_domain_);
                for (size_t a = 0; a < this->codeword_domain_size_; ++a)
                {
                    evals_of_response_poly[a] += current_evaluations[a] *
                                                 (input_matrix[i][a] + addition_matrix[i][a]+target_matrix[i][a]);// s()*f_w()
                }
                for(std::size_t k=0;k<this->num_queries_;k++){
                    query_s_matrix[i][k]=current_evaluations[query_set[k]];
                }
            }

            // 构造内积
            for(std::size_t k=0;k<this->num_queries_;k++){
                std::size_t pos=query_set[k];
                std::vector<FieldT> secret_vec;
                std::vector<FieldT> public_vec;
                secret_vec.resize(this->RS_col_systematic_domain_size_);
                public_vec.resize(this->RS_col_systematic_domain_size_);
                for(std::size_t i=0;i<this->RS_col_systematic_domain_size_;i++){
                    secret_vec[i]=input_matrix[i][pos]+addition_matrix[i][pos]+target_matrix[i][pos];
                    public_vec[i]=query_s_matrix[i][k];
                }
                polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                polynomial<FieldT> public_poly(IFFT_over_field_subset(public_vec,this->RS_col_systematic_domain_));
                secret_polys[count]=std::move(secret_poly);
                public_polys[count]=std::move(public_poly);
                target_sum+=evals_of_response_poly[pos];
                count++;
            }
        }

        const std::vector<FieldT> response_poly_coefficients =
                IFFT_over_field_subset<FieldT>(evals_of_response_poly, this->codeword_domain_);
        std::vector<FieldT> reduced_coefficients(&response_poly_coefficients[0],
                                                 &response_poly_coefficients[this->response_size_]);
        this->response_polys_coefficients[h]=reduced_coefficients;
    }


}


template<typename FieldT>
bool interleaved_lincheck_ot_protocol<FieldT>::verifier_predicate(const std::vector<FieldT> &supplementary_input,//v_bar
                                                                  std::vector<std::vector<FieldT>> random_linear_combinations,
                                                                  const std::vector<std::vector<FieldT>> &response_mes,
                                                                  const std::vector<std::size_t>& query_set)
{
    //verifier attaches random combination
    assert(!random_linear_combinations.empty());
    const std::vector<FieldT> codeword_elements = this->codeword_domain_.all_elements(); // eta

    for (size_t h = 0; h < this->num_interactions_; ++h) {
        /* EQUALITY TEST: does the polynomial that was sent sum to 0 over the systematic domain, as
           it should if the claimed statement is true? */

        libff::enter_block("Lincheck: equality test (p_0 sums to 0 over systematic domain)");
        std::vector<FieldT> response = response_mes[h];

        const std::vector<FieldT> evaluations = FFT_over_field_subset<FieldT>(response,
                                                                              this->extended_systematic_domain_);
        polynomial<FieldT> response_poly(std::move(response));

        FieldT equality_lhs(0);
        for (size_t d = 0; d < this->systematic_domain_size_; ++d) {
            const std::size_t idx = this->extended_systematic_domain_.reindex_by_subset(
                    this->systematic_domain_.dimension(), d);
            equality_lhs += evaluations[idx];
        }

        if (equality_lhs != FieldT(0)) {
            equality_lhs.print();
            return false;
        }
        libff::leave_block("Lincheck: equality test (p_0 sums to 0 over systematic domain)");
    }
    return true;
}

} // namespace ligero
