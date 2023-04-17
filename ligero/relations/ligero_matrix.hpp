#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "ligero/relations/variable.hpp"
namespace ligero{

template<typename FieldT>
using naive_sparse_matrix = std::vector<std::map<std::size_t, FieldT>>;

template<typename FieldT>
class ligero_multi_constraint{
public:

    linear_combination<FieldT> x_, y_, z_;

    ligero_multi_constraint() {};
    ligero_multi_constraint(const linear_combination<FieldT> &x,
                    const linear_combination<FieldT> &y,
                    const linear_combination<FieldT> &z);
};

template<typename FieldT>
class ligero_constraint_system{
public:
    std::size_t primary_input_size_;
    std::size_t auxiliary_input_size_;

    std::vector<ligero_multi_constraint<FieldT> > multi_constraints_;
    std::vector<linear_combination<FieldT>> add_constraints_;

    ligero_constraint_system() : primary_input_size_(0), auxiliary_input_size_(0) {}

    std::size_t num_inputs() const;
    std::size_t num_variables() const;
    std::size_t num_multi_constraints() const;
    std::size_t num_add_constraints() const;

    bool is_satisfied(const std::vector<FieldT> &primary_input,
                      const std::vector<FieldT> &auxiliary_input) const;
    bool is_satisfied(const std::vector<FieldT> &full_variable_assignment) const;

    void add_multi_constraint(const ligero_multi_constraint<FieldT> &c);
    void add_add_constraint(const linear_combination<FieldT> &c);

    naive_sparse_matrix<FieldT> X_matrix() const;
    naive_sparse_matrix<FieldT> Y_matrix() const;
    naive_sparse_matrix<FieldT> Z_matrix() const;
    naive_sparse_matrix<FieldT> ADD_matrix() const;

};
template<typename FieldT>
std::vector<FieldT> assignment_from_inputs(const std::vector<FieldT> &primary_input,
                                               const std::vector<FieldT> &auxiliary_input);

}//namespace

#include <ligero/relations/ligero_matrix.tcc>

#endif
