#ifndef Test_RScode_
#define Test_RScode_

#include <cstddef>
#include <memory>
#include <vector>

#include "ligero/algebra/polynomials/polynomial.hpp"
#include "ligero/algebra/field_subset/subspace.hpp"
#include "ligero/iop/iop.hpp"

namespace ligero{
template<typename FieldT>
class interleaved_RScode{
protected:
    iop_protocol<FieldT> &IOP_;
    domain_handle codeword_domain_handle_;
    domain_handle systematic_domain_handle_;
    domain_handle RS_col_systematic_domain_handle_;
    field_subset<FieldT> codeword_domain_;
    field_subset<FieldT> systematic_domain_;
    field_subset<FieldT> RS_col_systematic_domain_;

    std::size_t codeword_domain_size_;
    std::size_t systematic_domain_size_;
    std::size_t RS_col_systematic_domain_size_;
    std::size_t response_size_;
    std::size_t num_oracles_;
    std::size_t num_queries_;
    std::size_t num_col_queries_;
    std::size_t num_interactions_;

    bool make_zk_;
    field_subset_type field_subset_type_;

    std::vector<std::vector<std::vector<FieldT>>> query_col;// 打开的列
//    std::vector<std::size_t> query_pos;//打开的列号

    std::vector<verifier_random_message_handle> random_linear_combination_handles_;
    std::vector<prover_message_handle> response_handles_;

    std::vector<oracle_handle_ptr> blinding_vector_row_oracle_handles_;
    std::vector<oracle_handle_ptr> input_vector_row_oracle_handles_;
    std::vector<random_query_position_handle> query_position_handles_;
    std::vector<std::vector<query_handle>> input_vector_row_oracles_by_position_handles_;
    std::vector<std::vector<query_handle>> blinding_vector_row_oracles_by_position_handles_;
public:
    interleaved_RScode(iop_protocol<FieldT> &IOP,
                       const domain_handle &codeword_domain,//L
                       const domain_handle &systematic_domain,//H
                       const domain_handle &RS_col_systematic_domain_handle,
                       const std::size_t num_oracles, // 编码前矩阵的高度 是RS_sys_domain_size
                       const std::size_t num_queries,//打开的列
                       const std::size_t num_interactions,//进行的轮数
                       const bool make_zk,
                       const field_subset_type domain_type);
    void calculate_and_submit_responses(const std::vector<std::vector<FieldT>> &IRS_matrix,
                                        const std::vector<std::size_t>& query_set,
                                        const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                        const std::vector<std::vector<FieldT>> &blinding_matrix);

    void calculate_and_submit_responses(const std::vector<std::vector<FieldT>> &IRS_matrix,
                                        const std::vector<std::size_t>& query_set,
                                        const std::vector<std::vector<FieldT>> &random_linear_combinations,
                                        const std::vector<polynomial<FieldT>> public_polys,
                                        const std::vector<std::vector<FieldT>> &blinding_matrix);

    bool verifier_predicate(const std::vector<std::vector<FieldT>> &random_linear_combinations,
                            const std::vector<std::vector<FieldT>> &response_mes,
                            const std::vector<std::size_t> &query_pos);
    std::vector<std::vector<FieldT>> evals_of_response_polys;
    std::vector<polynomial<FieldT>> secret_polys;
    std::vector<polynomial<FieldT>> public_polys;
    FieldT target_sum;
};
}
#include "ligero/protocols/encoded/ligero/interactive_RS_test.tcc"

#endif