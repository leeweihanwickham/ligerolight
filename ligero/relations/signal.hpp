/**
****************************************************************************

 Description of iostream of a two-input arithmetic circuit gate

 a signal can be described as a_i or b_i
 a stands for a primary input
 b stands for an output of a gate

 *****************************************************************************/
#ifndef SIGNAL_HPP_
#define SIGNAL_HPP_

#include <cstddef>
namespace ligero{
typedef size_t s_index;// size_t==s_index == unsigned int
// 定义电路的线值
template<typename FieldT>
class signal{
    public:

        s_index index;
	bool is_primary_input;
	FieldT value;

	signal(const s_index index = 0,const bool b = true,FieldT t=FieldT(0)):index(index),is_primary_input(b),value(t){};

	void print();
};

}
#include "signal.tcc"

#endif
