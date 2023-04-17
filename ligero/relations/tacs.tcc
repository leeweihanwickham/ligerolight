#ifndef TACS_TCC_
#define TACS_TCC_

#include <vector>
#include <cassert>
#include<cstdio>
#include <queue>
#include <random>
namespace ligero{
template<typename FieldT>
bool tacs_gate<FieldT>::is_valid() const
{
    return(this->output.is_primary_input==false);
}

template<typename FieldT>
FieldT tacs_gate<FieldT>::evaluate(FieldT l_value,FieldT r_value)
{
    if(this->is_multi) {
        return l_value*r_value;
    }
    else return l_value+r_value;
}

template<typename FieldT>
std::size_t tacs_circuit<FieldT>::num_gates() const
{
    return gates.size();
}

template<typename FieldT>
std::size_t tacs_circuit<FieldT>::num_wires() const
{
    return circuit_input_size + num_gates();
}

//真正计算电路的输出
template<typename FieldT>
std::vector<FieldT> tacs_circuit<FieldT>::get_all_wires(const std::vector<FieldT> &input,const std::size_t RealInputNum){
    assert(input.size()==circuit_input_size);
    std::vector<FieldT> all_wires;
//    std::vector<treepoint<FieldT>> treeQ;
//    int flag;
    // input 是初始电路输入 2^k 按输入顺序排列a1 a2 a3 a4
    // 向量的插入  在all_wires(一开始是空的)的末尾插入input向量  input就是电路的输入
    all_wires.insert(all_wires.end(),input.begin(),input.end());

    for(std::size_t i=0,j=0;i<num_gates();i++){
//       treepoint<FieldT> treeP;
        FieldT l,r;

// 这里的index是在example_tacs.tcc里定义
// index从0开始 分别计算输入索引(a_i) 和 中间值索引(b_i)

        if(gates[i].l_input.is_primary_input){
            l = all_wires[gates[i].l_input.index];
        }
        else {l = all_wires[gates[i].l_input.index+circuit_input_size];}
        if(gates[i].r_input.is_primary_input)
            r = all_wires[gates[i].r_input.index];
        else r = all_wires[gates[i].r_input.index+circuit_input_size];

        FieldT gate_output = gates[i].evaluate(l,r);
        all_wires.emplace_back(gate_output);
//        std::cout<<gates[i].l_input.index<<" "<<gates[i].r_input.index<<" "<<gates[i].output.index<<" "<<l<<" "<<r<<" "<<gate_output<<std::endl;
    }

//  circuit_input_size=2^n  满二叉树的电路输入
    std::size_t NeedDeleteWiresNum;
    NeedDeleteWiresNum=circuit_input_size-RealInputNum;
    std::vector<FieldT> Real_wires;
    Real_wires.insert(Real_wires.end(),all_wires.begin()+NeedDeleteWiresNum*2,all_wires.end());
    std::vector<tacs_gate<FieldT>> Realgates;
    Realgates.insert(Realgates.end(),gates.begin()+NeedDeleteWiresNum,gates.end());


//  下面开始修改电路的索引等各种变量
    circuit_input_size=RealInputNum;
    all_wires=Real_wires;
    gates=Realgates;
//    修改gate里signal索引  输入门个数 circuit_input_size/2 - NeedDeleteNum/2  == RealInputNum/2
//    首先处理输入门  判断奇数个输入还是偶数个输入  奇数则说明最后一个输入门左输入是初始输入  右输入不是
    if((RealInputNum&0x1)==1){
        //奇数个输入  一共RealInputNum/2+1个
        std::size_t i,wirecounter=0,j=0;
        for(i=0;i<RealInputNum/2;i++){
            gates[i].output.index=i;
            gates[i].l_input.index=wirecounter;
            gates[i].r_input.index=wirecounter+1;
            gates[i].l_input.is_primary_input=true;
            gates[i].r_input.is_primary_input= true;
            wirecounter+=2;
        }
        gates[i].output.index=i;
        gates[i].l_input.is_primary_input=true;
        gates[i].r_input.is_primary_input= false;
        gates[i].l_input.index=gates[i-1].r_input.index+1;
        gates[i].r_input.index=j;
        i++;
        j++;
//        处理非输入门
        for(;i<gates.size();i++){
           gates[i].l_input.index=j;
           gates[i].r_input.index=j+1;
           gates[i].output.index=i;
           j+=2;
        }
    } else{
        std::size_t i,j=0,k=0;
        for(i=0;i<RealInputNum/2;i++){
            gates[i].output.index=i;
            gates[i].l_input.index=k;
            gates[i].r_input.index=k+1;
            gates[i].l_input.is_primary_input=true;
            gates[i].r_input.is_primary_input=true;
            k+=2;
        }
        for(;i<gates.size();i++){
            gates[i].output.index=i;
            gates[i].l_input.index=j;
            gates[i].r_input.index=j+1;
            j+=2;
        }
    }

    assert(all_wires.size()==num_wires());
    return all_wires;
}

template<typename FieldT>
std::vector<FieldT> tacs_circuit<FieldT>::get_random_wires(const std::vector<FieldT> &input)
{
    return circuit_input_size + num_gates();
}


template<typename FieldT>
bool tacs_circuit<FieldT>::is_legal() const
{
    for(std::size_t i=0;i<num_gates()-1;i=i+1){
	if(gates[i].is_circuit_output == true)return false;
    }
    if(gates[num_gates()-1].is_circuit_output == false) return false;
    else return true;
}

template<typename FieldT>
void tacs_circuit<FieldT>::add_gate(tacs_gate<FieldT> &g){
    assert(g.output.index == num_gates());
    assert(!((g.l_input.is_primary_input == false)&&(g.l_input.index>=num_gates())));
    assert(!((g.r_input.is_primary_input == false)&&(g.r_input.index>=num_gates())));
    assert(g.is_vaild());
    gates.emplace_back(g);
    if(g.l_input.is_primary_input){
	circuit_input_size++;
    }
    if(g.r_input.is_primary_input){
	circuit_input_size++;
    }
}

template<typename FieldT>
void tacs_circuit<FieldT>::generate_gate(signal<FieldT> l,signal<FieldT> r,bool type){
    signal<FieldT> o(num_gates(),false);
    assert(!((l.is_primary_input == false)&&(l.index>=num_gates())));
    assert(!((r.is_primary_input == false)&&(r.index>=num_gates())));
    if(l.is_primary_input){
	circuit_input_size++;
    }
    if(r.is_primary_input){
	circuit_input_size++;
    }
    tacs_gate<FieldT> g(l,r,o,false,type);
    // 在gates末尾添加门g 输出为o index是num_gates.size()
    gates.emplace_back(g);
}

template<typename FieldT>
void tacs_circuit<FieldT>::circuit_encapsulation(){
    gates[num_gates()-1].is_circuit_output = true;
}

template<typename FieldT>
void tacs_circuit<FieldT>::print(){
    for(std::size_t i=0;i<num_gates();i=i+1){
        gates[i].l_input.print();
        gates[i].r_input.print();
        gates[i].output.print();
        if(gates[i].is_multi)printf("*");
        else printf("+");
        printf("\n");
    }
}

}
#endif
