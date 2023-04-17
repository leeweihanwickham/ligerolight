#ifndef TACS_HPP_
#define TACS_HPP_

#include <cstddef>
#include <vector>
#include "signal.hpp"

/*Two-inputs Arithmetic Circuit Satisfaction(TACS)
*
*A TACS gate can calculate (a_1 + a_2) or (a_1 * a_2) on a certain field.
*
*A TACS circuit describes a map C:|F|^n->|F|.n represents size of primary input of the circuit.
*
*/
namespace ligero{
template<typename FieldT>
struct tacs_gate{
    signal<FieldT> l_input;
    signal<FieldT> r_input;
    signal<FieldT> output;

    bool is_circuit_output;
    bool is_multi;

    tacs_gate(signal<FieldT> l,signal<FieldT> r,signal<FieldT> o,bool b_1 = false,bool b_2 = false){
        l_input = l;
        r_input = r;
        output = o;
        is_circuit_output = b_1;
        is_multi = b_2;
    }

    bool is_valid() const;
    // 就是域的计算 加法或乘法
    FieldT evaluate(FieldT l_value,FieldT r_value);
};
//template<typename FieldT>
//struct treepoint{
//    FieldT lnum1,rnum1,outnum1;
//    struct treepoint *outpoint1;
//    struct treepoint *lpoint1;
//    struct treepoint *rpiont1;
//};
template<typename FieldT>
class tacs_circuit{
    public:
        std::size_t circuit_input_size;
        std::vector<tacs_gate<FieldT>> gates;

        tacs_circuit():circuit_input_size(0){};

        std::size_t num_gates() const;
        std::size_t num_wires() const;//num_gates()+num_inputs()
	
	std::vector<FieldT> get_all_wires(const std::vector<FieldT> &input,std::size_t Realnum);
	bool is_legal() const;

	void add_gate(tacs_gate<FieldT> &g);
	void generate_gate(signal<FieldT> l,signal<FieldT> r,bool type);
    std::vector<FieldT> get_random_wires(const std::vector<FieldT> &input);

	void circuit_encapsulation();

	void print();
};

}//namespace
#include "tacs.tcc"

#endif
