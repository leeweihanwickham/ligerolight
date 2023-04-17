#ifndef EXAMPLE_TACS_HPP_
#define EXAMPLE_TACS_HPP_

#include "ligero/relations/tacs.hpp"

namespace ligero{

template<typename FieldT>
struct tacs_example{
    tacs_circuit<FieldT> random_circuit;
    std::vector<FieldT> primary_input;
    std::vector<FieldT> auxiliary_input;

    tacs_example<FieldT>(const tacs_circuit<FieldT> &rc,
                         const std::vector<FieldT> &primary,
                         const std::vector<FieldT> &auxiliary):
	random_circuit(rc),
	primary_input(primary),
	auxiliary_input(auxiliary){};
    tacs_example<FieldT>(const tacs_circuit<FieldT> &&rc,
                         const std::vector<FieldT> &&primary,
                         const std::vector<FieldT> &&auxiliary):
	random_circuit(std::move(rc)),
	primary_input(std::move(primary)),
	auxiliary_input(std::move(auxiliary)){};
};
template<typename FieldT>
tacs_example<FieldT> generate_tacs_example(const int dim,std::size_t RealInputnum,std::size_t MultiGatesNum);
//Generate a arithmetic circuit which likes a full binary tree.Primary input number equals 2^dim
}
#include "example_tacs.tcc"

#endif
