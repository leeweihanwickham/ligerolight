#ifndef MATRIX_GEN_TCC_
#define MATRIX_GEN_TCC_

#include "tacs.hpp"
#include "signal.hpp"
#include "variable.hpp"

namespace ligero{

template<typename FieldT>
ligero_constraint_system<FieldT> ligero_constraint_gen(const tacs_circuit<FieldT> &t){
    ligero_constraint_system<FieldT> ls;

    ls.primary_input_size_ = t.circuit_input_size;
    ls.auxiliary_input_size_ = t.gates.size();

    for(int i=0;i<t.num_gates();i++){
        linear_combination<FieldT> X,Y,Z,ADD;
        if(t.gates[i].is_multi){
            variable<FieldT> vx,vy;
            if(t.gates[i].l_input.is_primary_input){
                vx.index_ = t.gates[i].l_input.index+1;
            }
            else vx.index_ = t.gates[i].l_input.index+t.circuit_input_size+1;
            linear_term<FieldT> x(vx,FieldT(1));
            X.terms.emplace_back(x);

            if(t.gates[i].r_input.is_primary_input){
                vy.index_ = t.gates[i].r_input.index+1;
            }
            else vy.index_ =t.gates[i].r_input.index+t.circuit_input_size+1;
            linear_term<FieldT> y(vy,FieldT(1));
            Y.terms.emplace_back(y);

            variable<FieldT> vz(t.gates[i].output.index+t.circuit_input_size+1);
            linear_term<FieldT> z(vz,FieldT(1));
            Z.terms.emplace_back(z);

            ls.add_multi_constraint(ligero_multi_constraint<FieldT>(X,Y,Z));
        }
        else{
            variable<FieldT> vx,vy;
            if(t.gates[i].l_input.is_primary_input){
                vx.index_ = t.gates[i].l_input.index+1;
            }
            else vx.index_ = t.gates[i].l_input.index+t.circuit_input_size+1;
            linear_term<FieldT> x(vx,FieldT(1));
            ADD.terms.emplace_back(x);

            if(t.gates[i].r_input.is_primary_input){
                vy.index_ = t.gates[i].r_input.index+1;
            }
            else vy.index_ =t.gates[i].r_input.index+t.circuit_input_size+1;
            linear_term<FieldT> y(vy,FieldT(1));
            ADD.terms.emplace_back(y);

            variable<FieldT> vz(t.gates[i].output.index+t.circuit_input_size+1);
            linear_term<FieldT> z(vz,FieldT(1));
            ADD.terms.emplace_back(z);

            ls.add_add_constraint(ADD);
        }

    }
    return ls;
}

}//namespace

#endif
