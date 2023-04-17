#ifndef EXAMPLE_TACS_TCC_
#define EXAMPLE_TACS_TCC_

#include <random>
#include <queue>
#include "example_tacs.hpp"

namespace ligero{
template<typename FieldT>
tacs_example<FieldT> generate_tacs_example(const int dim,std::size_t realInputnum,std::size_t MultiGatesNum){
    tacs_circuit<FieldT> rc;
    std::queue<signal<FieldT>> q;
    std::size_t num_input = 1ull<<dim;  //2^dim
    std::vector<FieldT> primary_input;
//    std::cout<<"输入的个数"<<num_input<<std::endl;
    std::default_random_engine e;
    std::uniform_int_distribution<unsigned> u(0,1);//生成随机整数值 i ，均匀分布于闭区间 [0, 1]
//    产生输入值  队列q用来标识顺序  a是电路线值  true表示是初始输入
    for(std::size_t i=0;i<num_input;i++){
        // i作为sigal的index
        signal<FieldT> a(i,true);
        q.push(a);
        primary_input.emplace_back(FieldT::random_element());
        // 在primary_input尾部加入随机数 表示输入值   即primary_input才是真正记录输入的
    }

    std::vector<bool> Deletegatestype(num_input-realInputnum, false);//需要删掉的门都是加法门
    std::vector<bool> Multigatestype(MultiGatesNum,true);
    std::vector<bool> Addgatestype(realInputnum-1-MultiGatesNum, false);
    std::vector<bool> gatestype;
    gatestype.insert(gatestype.end(),Multigatestype.begin(),Multigatestype.end());
    gatestype.insert(gatestype.end(),Addgatestype.begin(),Addgatestype.end());
    std::shuffle(gatestype.begin(),gatestype.end(), std::mt19937(std::random_device()()));
    gatestype.insert(gatestype.begin(),Deletegatestype.begin(),Deletegatestype.end());

    std::size_t count = 0;
    while(q.size()>1){
        signal<FieldT> l,r;
        l = q.front();//获取队列的第一个元素 并不具有电路的值 而是一个标志 代表某个输入
        q.pop();
        r = q.front();
        q.pop();
        // generate_gate方法产生一个门 并将这个门推入到队列 gates 的末尾  因此gates队列的门是按顺序排列的
        // 通过l r的输入产生一个门 将输出标记为b_i 并推入到队列中
//        rc.generate_gate(l,r,(u(e)==1));
        rc.generate_gate(l,r,gatestype[count]);
        // 队列最前面的是输入向量 当输入向量遍历结束后 就产生了输入门 此时队列中是按照(b0 b1 b2 b3 ...)排列的 其中b0=a0+a1
        // 这里的output也是没有具体的值 只是一个signal
        q.push(rc.gates[count].output);
        count++;
    }


    // 真正计算电路的输出  事实上vector<gates>里就是signal 没有具体的值
    num_input=realInputnum;
    std::vector<FieldT> w = rc.get_all_wires(primary_input,num_input);
//    auxiliary_input  非输入的线值

    std::vector<FieldT> auxiliary_input(w.begin()+num_input,w.end());
    primary_input.clear();
    primary_input.insert(primary_input.end(),w.begin(),w.begin()+num_input);

    return tacs_example<FieldT>(std::move(rc),std::move(primary_input),std::move(auxiliary_input));
}

}
#endif
