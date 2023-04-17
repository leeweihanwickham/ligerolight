#include <stdexcept>
#include <thread>
#include <libff/common/profiling.hpp>
#include "ligero/algebra/fft.hpp"
namespace ligero {

template<typename FieldT>
interleaved_rowcheck_protocol<FieldT>::interleaved_rowcheck_protocol(
    iop_protocol<FieldT> &IOP,
    const domain_handle &codeword_domain_handle,
    const domain_handle &systematic_domain_handle,
    const domain_handle &extended_systematic_domain_handle,
    const domain_handle &RS_col_systematic_domain_handle,
    const std::size_t num_oracles_input,
    const std::size_t num_oracles,
    const std::size_t num_queries,
    const std::size_t num_interactions,
    const bool make_zk,
    const field_subset_type domain_type) :
    IOP_(IOP),
    codeword_domain_handle_(codeword_domain_handle),
    systematic_domain_handle_(systematic_domain_handle),
    RS_col_systematic_domain_handle_(RS_col_systematic_domain_handle),
    num_oracles_input_(num_oracles_input),
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

    this->extended_systematic_domain_ = this->IOP_.get_domain(extended_systematic_domain_handle);

    this->RS_col_systematic_domain_=this->IOP_.get_domain(RS_col_systematic_domain_handle);
    this->RS_col_systematic_domain_size_=this->RS_col_systematic_domain_.num_elements();
    this->response_size_ = 2 * this->systematic_domain_size_;
}


/* Proving */
template<typename FieldT>
void interleaved_rowcheck_protocol<FieldT>::calculate_and_submit_responses(const std::vector<std::vector<FieldT>> &Ux_matrix,
                                                                           const std::vector<std::vector<FieldT>> &Uy_matrix,
                                                                           const std::vector<std::vector<FieldT>> &Uz_matrix,
                                                                           const std::vector<std::size_t>& query_set,
                                                                           const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                                                           const std::vector<std::vector<FieldT>> &blinding_vectors)
{

//    std::vector<std::vector<FieldT>> response_coefficients;
    assert(!random_linear_combinations.empty());
    this->response_coefficients.resize(this->num_interactions_);
    target_sum=FieldT(0);
    std::size_t count=0;
    secret_polys.resize(this->num_interactions_*this->num_queries_);
    public_polys.resize(this->num_interactions_*this->num_queries_);
    for (size_t h = 0; h < this->num_interactions_; ++h)
    {
        std::vector<FieldT> evals_of_response_poly(this->codeword_domain_size_, FieldT(0)); // evaluations of p_0

        /** Build the response polynomial's evaluations row by row. */
        // 考虑到Uw的blinding位置在(m2)  Uxyz的高度是RS_col_sys_size
        if(this->make_zk_){
            for (std::size_t i = this->num_oracles_input_+1; i < this->num_oracles_input_+1+this->num_oracles_; ++i)
            {
                /** For each column, add to that column's corresponding response polynomial evaluation
                 *      r * (p_x[col] * p_y[col] - p_z[col])
                 *  from this row. */
                for (size_t j = 0; j < this->codeword_domain_size_; ++j)
                {
                    FieldT val = (Ux_matrix[i][j] * Uy_matrix[i][j])
                                 - Uz_matrix[i][j];
                    evals_of_response_poly[j] += random_linear_combinations[h][i] * val;
                }
            }
            // blinding
            for (std::size_t j = 0; j < this->codeword_domain_size_; ++j)
            {
                evals_of_response_poly[j] += random_linear_combinations[h][this->num_oracles_input_+1+this->num_oracles_] * blinding_vectors[h][j];
            }

            //构造内积
            for(unsigned long pos : query_set){
                std::vector<FieldT> secret_vec;
                secret_vec.resize(this->RS_col_systematic_domain_size_);
                for(std::size_t i=0;i<this->RS_col_systematic_domain_size_;i++){
                    secret_vec[i]=Ux_matrix[i][pos]*Uy_matrix[i][pos]-Uz_matrix[i][pos];
                }
                secret_vec[this->num_oracles_input_+1+this->num_oracles_]=blinding_vectors[h][pos];

                if (h==0){
                    polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                    polynomial<FieldT> public_poly(IFFT_over_field_subset(random_linear_combinations[h],this->RS_col_systematic_domain_));
                    secret_polys[count]=std::move(secret_poly);
                    public_polys[count]=std::move(public_poly);
                    target_sum+=evals_of_response_poly[pos];
                    count++;
                }
                else{
                    polynomial<FieldT> public_poly(IFFT_over_field_subset(random_linear_combinations[h],this->RS_col_systematic_domain_));
                    secret_polys[count]=secret_polys[count-this->num_queries_];
                    public_polys[count]=std::move(public_poly);
                    target_sum+=evals_of_response_poly[pos];
                    count++;
                }

            }
        } else{
            for (size_t i = 0; i < this->RS_col_systematic_domain_size_; ++i)
            {
                /** For each column, add to that column's corresponding response polynomial evaluation
                 *      r * (p_x[col] * p_y[col] - p_z[col])
                 *  from this row. */
                for (size_t j = 0; j < this->codeword_domain_size_; ++j)
                {
                    FieldT val = (Ux_matrix[i][j] * Uy_matrix[i][j])
                                 - Uz_matrix[i][j];
                    evals_of_response_poly[j] += random_linear_combinations[h][i] * val;
                }
            }
            //构造内积
            for(unsigned long pos : query_set){
                std::vector<FieldT> secret_vec;
                secret_vec.resize(this->RS_col_systematic_domain_size_);
                for(std::size_t i=0;i<this->RS_col_systematic_domain_size_;i++){
                    secret_vec[i]=Ux_matrix[i][pos]*Uy_matrix[i][pos]-Uz_matrix[i][pos];
                }
                polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                polynomial<FieldT> public_poly(IFFT_over_field_subset(random_linear_combinations[h],this->RS_col_systematic_domain_));
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
        this->response_coefficients[h]=reduced_coefficients;
    }

}

template<typename FieldT>
void interleaved_rowcheck_protocol<FieldT>::calculate_and_submit_responses(const std::vector<std::vector<FieldT>> &Ux_matrix,
                                                                           const std::vector<std::vector<FieldT>> &Uy_matrix,
                                                                           const std::vector<std::vector<FieldT>> &Uz_matrix,
                                                                           const std::vector<std::size_t>& query_set,
                                                                           const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                                                           const std::vector<polynomial<FieldT>> input_public_polys,
                                                                           const std::vector<std::vector<FieldT>> &blinding_vectors)
{

//    std::vector<std::vector<FieldT>> response_coefficients;
    assert(!random_linear_combinations.empty());
    this->response_coefficients.resize(this->num_interactions_);
    target_sum=FieldT(0);
    std::size_t count=0;
    secret_polys.resize(this->num_interactions_*this->num_queries_);
    public_polys.resize(this->num_interactions_*this->num_queries_);
    for (size_t h = 0; h < this->num_interactions_; ++h)
    {
        std::vector<FieldT> evals_of_response_poly(this->codeword_domain_size_, FieldT(0)); // evaluations of p_0

        /** Build the response polynomial's evaluations row by row. */
        // 考虑到Uw的blinding位置在(m2)  Uxyz的高度是RS_col_sys_size
        if(this->make_zk_){
            for (std::size_t i = this->num_oracles_input_+1; i < this->num_oracles_input_+1+this->num_oracles_; ++i)
            {
                /** For each column, add to that column's corresponding response polynomial evaluation
                 *      r * (p_x[col] * p_y[col] - p_z[col])
                 *  from this row. */
                for (size_t j = 0; j < this->codeword_domain_size_; ++j)
                {
                    FieldT val = (Ux_matrix[i][j] * Uy_matrix[i][j])
                                 - Uz_matrix[i][j];
                    evals_of_response_poly[j] += random_linear_combinations[h][i] * val;
                }
            }
            // blinding
            for (std::size_t j = 0; j < this->codeword_domain_size_; ++j)
            {
                evals_of_response_poly[j] += random_linear_combinations[h][this->num_oracles_input_+1+this->num_oracles_] * blinding_vectors[h][j];
            }

            //构造内积
            for(unsigned long pos : query_set){
                if (h==0){
                    std::vector<FieldT> secret_vec;
                    secret_vec.resize(this->RS_col_systematic_domain_size_);
                    for(std::size_t i=0;i<this->RS_col_systematic_domain_size_;i++){
                        secret_vec[i]=Ux_matrix[i][pos]*Uy_matrix[i][pos]-Uz_matrix[i][pos];
                    }
                    secret_vec[this->num_oracles_input_+1+this->num_oracles_]=blinding_vectors[h][pos];
                    polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                    secret_polys[count]=std::move(secret_poly);
                    public_polys[count]=std::move(input_public_polys[count]);
                    target_sum+=evals_of_response_poly[pos];
                    count++;
                }
                else{
                    polynomial<FieldT> public_poly(IFFT_over_field_subset(random_linear_combinations[h],this->RS_col_systematic_domain_));
                    secret_polys[count]=secret_polys[count-this->num_queries_];
                    public_polys[count]=std::move(input_public_polys[count]);
                    target_sum+=evals_of_response_poly[pos];
                    count++;
                }

            }
        } else{
            for (size_t i = 0; i < this->RS_col_systematic_domain_size_; ++i)
            {
                /** For each column, add to that column's corresponding response polynomial evaluation
                 *      r * (p_x[col] * p_y[col] - p_z[col])
                 *  from this row. */
                for (size_t j = 0; j < this->codeword_domain_size_; ++j)
                {
                    FieldT val = (Ux_matrix[i][j] * Uy_matrix[i][j])
                                 - Uz_matrix[i][j];
                    evals_of_response_poly[j] += random_linear_combinations[h][i] * val;
                }
            }
            //构造内积
            for(unsigned long pos : query_set){
                std::vector<FieldT> secret_vec;
                secret_vec.resize(this->RS_col_systematic_domain_size_);
                for(std::size_t i=0;i<this->RS_col_systematic_domain_size_;i++){
                    secret_vec[i]=Ux_matrix[i][pos]*Uy_matrix[i][pos]-Uz_matrix[i][pos];
                }
                polynomial<FieldT> secret_poly(IFFT_over_field_subset(secret_vec,this->RS_col_systematic_domain_));
                polynomial<FieldT> public_poly(IFFT_over_field_subset(random_linear_combinations[h],this->RS_col_systematic_domain_));
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
        this->response_coefficients[h]=reduced_coefficients;
    }

}

/* Verification */
template<typename FieldT>
bool interleaved_rowcheck_protocol<FieldT>::verifier_predicate(const std::vector<std::vector<FieldT>> &random_vectors,
                                                               const std::vector<std::vector<FieldT>> &response_mes,
                                                               std::vector<std::size_t> query_set)
{
    const std::vector<FieldT> codeword_elements = this->codeword_domain_.all_elements(); // eta

    for (size_t h = 0; h < this->num_interactions_; ++h)
    {
        /* EQUALITY TEST: does the polynomial that was sent equal 0 over the entire systematic
           domain, as it should if the claimed statement is true? */

        libff::enter_block("Rowcheck: equality test (p_0 equals 0 within systematic domain)");
        std::vector<FieldT> response = response_mes[h];
        const std::vector<FieldT> evaluations = FFT_over_field_subset<FieldT>(response, this->extended_systematic_domain_);
        const polynomial<FieldT> response_poly(std::move(response));

        // TODO zk
        for (size_t d = 0; d < this->systematic_domain_size_ ; ++d)
        {
            const std::size_t idx = this->extended_systematic_domain_.reindex_by_subset(
                    this->systematic_domain_.dimension(), d);
            if (evaluations[idx] != FieldT(0))
            {
                return false;
            }
        }
        // TODO zk
        libff::leave_block("Rowcheck: equality test (p_0 equals 0 within systematic domain)");
    }

    return true;
}

} // namespace ligero

