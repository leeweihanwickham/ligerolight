#ifndef TACS_TO_R1CS_HPP_
#define TACS_TO_R1CS_HPP_

#include <map>
#include <unordered_map>
#include <vector>
#include "tacs.hpp"
#include "r1cs.hpp"

namespace ligero{

std::unordered_map<std::size_t,std::size_t> umap;

template<typename FieldT>
void circuit_check(tacs_circuit<FieldT> &t,std::vector<FieldT> &p_i,std::vector<FieldT> &a_i);

//function "check"
template<typename FieldT>
void check(std::vector<std::size_t> &s,const tacs_circuit<FieldT> &t,const signal<FieldT> &v);

template<typename FieldT>
r1cs_constraint_system<FieldT> tacs_to_r1cs(const tacs_circuit<FieldT> &t);

}
#include "tacs_to_r1cs.tcc"

#endif
