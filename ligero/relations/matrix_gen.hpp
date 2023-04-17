#ifndef MATRIX_GEN_HPP_
#define MATRIX_GEN_HPP_

#include "tacs.hpp"
#include "ligero/relations/ligero_matrix.hpp"

namespace ligero{

template<typename FieldT>
ligero_constraint_system<FieldT> ligero_constraint_gen(const tacs_circuit<FieldT> &t);

}//namespace

#include <ligero/relations/matrix_gen.tcc>

#endif
