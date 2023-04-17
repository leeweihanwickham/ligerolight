#ifndef MATRIX_TCC_
#define MATRIX_TCC_

#include <algorithm>
#include <cassert>
#include <set>

#include <libff/common/utils.hpp>

namespace ligero {

template<typename FieldT>
ligero_multi_constraint<FieldT>::ligero_multi_constraint(const linear_combination<FieldT> &x,
                                         const linear_combination<FieldT> &y,
                                         const linear_combination<FieldT> &z) :
    x_(x), y_(y), z_(z)
{}
template<typename FieldT>
std::size_t ligero_constraint_system<FieldT>::num_inputs() const
{
    return primary_input_size_;
}

template<typename FieldT>
std::size_t ligero_constraint_system<FieldT>::num_variables() const
{
    return primary_input_size_ + auxiliary_input_size_;
}


template<typename FieldT>
std::size_t ligero_constraint_system<FieldT>::num_multi_constraints() const
{
    return multi_constraints_.size();
}
template<typename FieldT>
std::size_t ligero_constraint_system<FieldT>::num_add_constraints() const
{
    return add_constraints_.size();
}

template<typename FieldT>
bool ligero_constraint_system<FieldT>::is_satisfied(const std::vector<FieldT> &primary_input,
                                                    const std::vector<FieldT> &auxiliary_input) const
{
    assert(primary_input.size() == num_inputs());
    assert(primary_input.size() + auxiliary_input.size() == num_variables());

    const std::vector<FieldT> full_variable_assignment =
        assignment_from_inputs(primary_input, auxiliary_input);
    return this->is_satisfied(full_variable_assignment);
}

template<typename FieldT>
bool ligero_constraint_system<FieldT>::is_satisfied(const std::vector<FieldT> &full_variable_assignment) const
{
    for (size_t c = 0; c < multi_constraints_.size(); ++c)
    {
        const FieldT ares = multi_constraints_[c].x_.evaluate(full_variable_assignment);
        const FieldT bres = multi_constraints_[c].y_.evaluate(full_variable_assignment);
        const FieldT cres = multi_constraints_[c].z_.evaluate(full_variable_assignment);

        if (!(ares*bres == cres))
        {
            return false;
        }
    }
    for (size_t c = 0; c < add_constraints_.size(); ++c)
    {
        const FieldT zerores = add_constraints_[c].evaluate(full_variable_assignment);

        if (!(zerores == FieldT(0)))
        {
            return false;
        }
    }
    return true;
}

template<typename FieldT>
void ligero_constraint_system<FieldT>::add_multi_constraint(const ligero_multi_constraint<FieldT> &c)
{
    multi_constraints_.emplace_back(c);
}

template<typename FieldT>
naive_sparse_matrix<FieldT> ligero_constraint_system<FieldT>::X_matrix() const
{
    naive_sparse_matrix<FieldT> matrix;
    for (std::size_t i = 0; i < this->num_multi_constraints(); ++i)
    {
        const ligero_multi_constraint<FieldT> &constraint = this->multi_constraints_[i];

        std::vector<linear_term<FieldT>> terms = constraint.x_.terms;
        std::map<std::size_t, FieldT> values;
        for (std::size_t j = 0; j < terms.size(); ++j)
        {
            values.insert(std::pair<std::size_t, FieldT>(terms[j].index_, terms[j].coeff_));
        }

        matrix.push_back(std::move(values));
    }
    return matrix;
}
template<typename FieldT>
naive_sparse_matrix<FieldT> ligero_constraint_system<FieldT>::Y_matrix() const
{
    naive_sparse_matrix<FieldT> matrix;
    for (std::size_t i = 0; i < this->num_multi_constraints(); ++i)
    {
        const ligero_multi_constraint<FieldT> &constraint = this->multi_constraints_[i];

        std::vector<linear_term<FieldT>> terms = constraint.y_.terms;
        std::map<std::size_t, FieldT> values;
        for (std::size_t j = 0; j < terms.size(); ++j)
        {
            values.insert(std::pair<std::size_t, FieldT>(terms[j].index_, terms[j].coeff_));
        }

        matrix.push_back(std::move(values));
    }
    return matrix;
}
template<typename FieldT>
naive_sparse_matrix<FieldT> ligero_constraint_system<FieldT>::Z_matrix() const
{
    naive_sparse_matrix<FieldT> matrix;
    for (std::size_t i = 0; i < this->num_multi_constraints(); ++i)
    {
        const ligero_multi_constraint<FieldT> &constraint = this->multi_constraints_[i];

        std::vector<linear_term<FieldT>> terms = constraint.z_.terms;
        std::map<std::size_t, FieldT> values;
        for (std::size_t j = 0; j < terms.size(); ++j)
        {
            values.insert(std::pair<std::size_t, FieldT>(terms[j].index_, terms[j].coeff_));
        }

        matrix.push_back(std::move(values));
    }
    return matrix;
}
template<typename FieldT>
naive_sparse_matrix<FieldT> ligero_constraint_system<FieldT>::ADD_matrix() const
{
    naive_sparse_matrix<FieldT> matrix;
    for (std::size_t i = 0; i < this->num_add_constraints(); ++i)
    {

        std::vector<linear_term<FieldT>> terms = this->add_constraints_[i].terms;
        std::map<std::size_t, FieldT> values;
        for (std::size_t j = 0; j < terms.size(); ++j)
        {
            values.insert(std::pair<std::size_t, FieldT>(terms[j].index_, terms[j].coeff_));
        }

        matrix.push_back(std::move(values));
    }
    return matrix;
}

template<typename FieldT>
void ligero_constraint_system<FieldT>::add_add_constraint(const linear_combination<FieldT> &c)
{
    add_constraints_.emplace_back(c);
}

template<typename FieldT>
std::vector<FieldT> assignment_from_inputs(const std::vector<FieldT> &primary_input,
                                                                 const std::vector<FieldT> &auxiliary_input)
{
    std::vector<FieldT> full_variable_assignment = primary_input;
    full_variable_assignment.insert(full_variable_assignment.end(), auxiliary_input.begin(), auxiliary_input.end());
    return full_variable_assignment;
}

}//namespace

#endif
