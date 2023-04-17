#ifndef LIGERO_IOP_HPP_
#define LIGERO_IOP_HPP_

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>
#include <cmath>

#include "ligero/iop/iop.hpp"
#include "ligero/protocols/encoded/ligero/ligero.hpp"


#include "ligero/protocols/ldt/ldt_reducer.hpp"
#include "ligero/protocols/ldt/direct_ldt/direct_ldt.hpp"
#include "ligero/relations/tacs.hpp"
#include "ligero/relations/ligero_matrix.hpp"
#include "ligero/relations/matrix_gen.hpp"

namespace ligero{
template<typename FieldT>
class ligero_iop_parameters{
protected:
    std::size_t security_parameter_;
    std::size_t RS_extra_dimensions_;

    bool make_zk_;
    std::size_t zk_queries;
    float height_width_ratio_;
    field_subset_type domain_type_;

    std::size_t num_variables_;
    std::size_t num_constraints_;

    std::size_t num_oracles_input_;
    std::size_t num_oracle_vectors_;
    std::size_t systematic_domain_dim_;//log h
    std::size_t codeword_domain_dim_;//log l
    void set_soundness_parameters();
    void set_encoded_ligero_interactions(const size_t interactive_soundness_bits);
    void set_queries(const size_t query_soundness_bits);
    void configure_encoded_ligero_params(const size_t num_variables,
                                         const size_t num_constraints);

public:
    ligero_iop_parameters(const size_t security_parameter,
                          const size_t RS_extra_dimensions,
                          const size_t RS_col_extra_dimensions,
                          const std::vector<std::size_t> localization_parameter_array,
            // 描述了待编码前的谕示矩阵的横纵比
                          const float height_width_ratio,
                          const bool make_zk,
                          const field_subset_type domain_type,
            // 乘法门的约束数目，也就是x,y,z的长度
                          const size_t num_constraints,
            // 总变量数，也就是w的长度，P_add w = 0 , P_x w = x, P_y w = y, P_z w = z
                          const size_t num_variables,
                          const size_t parallel_opt);
    std::size_t systematic_domain_dim() const;
    std::size_t RS_extra_dimensions() const;
    std::size_t col_RS_extra_dimensions() const;
    bool make_zk() const;
    field_subset_type domain_type() const;
    void print() const;
    std::size_t parallel_opt_;
    encoded_ligero_parameters encoded_ligero_params_;
    std::size_t RS_col_extra_dimensions_;
    std::vector<std::size_t> localization_parameter_array_;
    std::size_t col_systematic_domain_dim_;
};


template<typename FieldT>
class ligero_iop {
protected:
    iop_protocol<FieldT> &IOP_;
    r1cs_constraint_system<FieldT> constraint_system_;
    ligero_iop_parameters<FieldT> parameters_;

    field_subset<FieldT> codeword_domain_;
    field_subset<FieldT> RS_col_codeword_domain_;
    std::shared_ptr<interleaved_r1cs_protocol<FieldT>> protocol_;

public:
    /* Initialization and registration */
    ligero_iop(iop_protocol<FieldT> &IOP,
                                const r1cs_constraint_system<FieldT> &constraint_system,
                                const ligero_iop_parameters<FieldT> &parameters);

    void register_interactions();
    void register_queries();

    /* Proving */
    float produce_oracle(const std::vector<FieldT> &primary_input,
                        const std::vector<FieldT> &auxiliary_input,
                        const std::map<std::size_t,FieldT> public_input);
    float produce_proof(const std::map<std::size_t,FieldT> public_input);
    /* Verification */
    bool verifier_predicate(const std::map<std::size_t,FieldT> public_input);
protected:
    void submit_random_blinding_vector(const oracle_handle_ptr &handle);
};
}
#include "ligero/protocols/ligero_iop.tcc"
#endif
