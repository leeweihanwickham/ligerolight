#ifndef TACS_TO_R1CS_TCC_
#define TACS_TO_R1CS_TCC_

#include "tacs.hpp"
#include "signal.hpp"
#include "variable.hpp"
namespace ligero{


template<typename FieldT>
// P_i  初始输入  a_i  非初始输入  t是电路
void circuit_check(tacs_circuit<FieldT> &t,std::vector<FieldT> &p_i,std::vector<FieldT> &a_i){
//    最后一个门如果不是乘法门  则增加了一个乘法门
    if(!t.gates[t.num_gates()-1].is_multi){
//      r:  is_primary_input-->true  作为右输入
        signal<FieldT> r(t.circuit_input_size,true);
        t.gates[t.num_gates()-1].is_circuit_output = false;
//        在gates里添加一个乘法门
        t.generate_gate(t.gates[t.num_gates()-1].output,r,true);
        t.circuit_encapsulation();
//        在初始输入p_i后面加一个值为1的输入  在a_i后面复制最后的输入
        p_i.emplace_back(FieldT(1));
        a_i.emplace_back(a_i.back());
    }
    std::vector<FieldT> r1cs_a_i;
    std::size_t pointer = 0;
//    把乘法门的结果取出来  这里的a_i和门序号是对应的
    for(std::size_t i=0;i<t.num_gates();i++){
        if(t.gates[i].is_multi){
            std::pair<std::size_t,std::size_t> temp(i+t.circuit_input_size,pointer+t.circuit_input_size);
            pointer++;
            umap.insert(temp);
            r1cs_a_i.emplace_back(a_i[i]);
        }
    }

    a_i = r1cs_a_i;

}

template<typename FieldT>
void check(std::vector<std::size_t> &s,const tacs_circuit<FieldT> &t,const signal<FieldT> &v){
    std::vector<signal<FieldT>> process;
    process.emplace_back(v);
    while(process.size()!=0){
        signal<FieldT> m=process.back();
        if(m.is_primary_input){
            s.emplace_back(m.index);
            process.pop_back();
        }
        //signal is a primary input
//        如果m不是初始输入  那么m.index代表门的index
        else{
            if(!t.gates[m.index].is_multi){
                process.pop_back();
                process.emplace_back(t.gates[m.index].l_input);
                process.emplace_back(t.gates[m.index].r_input);
            }
            else{
                s.emplace_back(m.index+t.circuit_input_size);
                process.pop_back();
            }
        }
    }
}

template<typename FieldT>
r1cs_constraint_system<FieldT> tacs_to_r1cs(const tacs_circuit<FieldT> &t){
    r1cs_constraint_system<FieldT> cs;
//    for(std::size_t i=0;i<t.gates.size();i++){
//        std::cout<<"l r o "<<i<<" "<<t.gates[i].l_input.index<<" "<<t.gates[i].r_input.index<<" "<<t.gates[i].output.index<<std::endl;
//    }
    std::vector<std::size_t> circuit_stack;
    std::vector<std::size_t> process;

    cs.primary_input_size_ = t.circuit_input_size;
    cs.auxiliary_input_size_ = umap.size();//乘法门的输入个数
    process.emplace_back(t.num_gates()-1);

    while(process.size()!=0){
        linear_combination<FieldT> A,B,C;
        std::size_t temp=process.back();
        std::size_t position = 0;
        check(circuit_stack,t,t.gates[temp].l_input);

        process.pop_back();
        for(std::size_t i=0;i<circuit_stack.size();i++){
            if(circuit_stack[i]>t.circuit_input_size-1){
                process.emplace_back(circuit_stack[i]-t.circuit_input_size);
                position = umap.at(circuit_stack[i]);
            }
            else position = circuit_stack[i];
            variable<FieldT> a(position+1);
            linear_term<FieldT> l(a,FieldT(1));
            A.terms.emplace_back(l);
        }
        circuit_stack.resize(0);
        check(circuit_stack,t,t.gates[temp].r_input);
        for(std::size_t i=0;i<circuit_stack.size();i++){
            if(circuit_stack[i]>t.circuit_input_size-1){
                process.emplace_back(circuit_stack[i]-t.circuit_input_size);
                position = umap.at(circuit_stack[i]);
            }
            else position = circuit_stack[i];
            variable<FieldT> a(position+1);
            linear_term<FieldT> l(a,FieldT(1));
            B.terms.emplace_back(l);
        }
        circuit_stack.resize(0);
        check(circuit_stack,t,t.gates[temp].output);
        for(std::size_t i=0;i<circuit_stack.size();i++){
            variable<FieldT> a(umap.at(circuit_stack[i])+1);
            linear_term<FieldT> l(a,FieldT(1));
            C.terms.emplace_back(l);
        }
        circuit_stack.resize(0);
        cs.add_constraint(r1cs_constraint<FieldT>(A, B, C));
    }
    return cs;
}

}

#endif
