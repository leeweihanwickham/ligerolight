#include <cstdint>
#include <stdexcept>
#include <ctime>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/fields/prime_base/fields_64.hpp>
#include "ligero/algebra/fft.hpp"
#include "ligero/algebra/polynomials/polynomial.hpp"
#include "ligero/algebra/field_subset/subspace.hpp"
#include <libff/common/utils.hpp>
#include "ligero/protocols/ligero_iop.hpp"
#include "ligero/relations/r1cs.hpp"
#include "ligero/relations/variable.hpp"
#include "ligero/bcs/merkle_tree.hpp"
#include "ligero/bcs/Newmerkle.hpp"

using namespace ligero;

std::vector<std::size_t> get_queries(std::size_t query_num,std::size_t domain_size){
    std::vector<std::size_t> query_set;
    for(std::size_t i=0;i<query_num;i++){
        bool is_repeat= true;
        std::size_t val;
        while(is_repeat){
            val=std::rand()%domain_size;

            std::vector<std::size_t>::iterator it;
            it= std::find(query_set.begin(),query_set.end(),val);
            if(it==query_set.end()){
                is_repeat= false;
            }
        }
        query_set.emplace_back(val);
    }
//    std::cout<<query_set.size()<<std::endl;
    for(unsigned long & i : query_set){
        i+=(domain_size-1);
//        std::cout<<i<<" ";
    }
//    std::cout<<"\n";
    return query_set;
}

int main(){
    typedef libff::Fields_64 FieldT;
    std::shared_ptr<merkle<FieldT>> merkleTree;
    const std::size_t width=1024;//矩阵宽度 2^k
    const std::size_t height=8;//矩阵高度 随意
    const bool type= true;// 按列做承诺     若按行做承诺 则将heigth改为2^k，width随意，type改为false
    merkleTreeParameter par,par2,par3;
    std::vector<std::vector<FieldT>> value_for_commit;
    value_for_commit.resize(height);

    for(std::size_t i=0;i<height;i++){
        value_for_commit[i]=random_FieldT_vector<FieldT>(width);
    }

    std::vector<std::size_t> pos=get_queries(4,width);

    // width is the
    merkleTree.reset(new merkle<FieldT>(
            width,
            pos,
            false
    ));

    // commit a matrix
    par=merkleTree->create_merklePar_of_matrix(value_for_commit);
    bool suc2=merkleTree->verify_merkle_commit(par);

    // commit a vector
    par2 = merkleTree->create_merklePar_of_vec(value_for_commit[1]);

    bool suc3 = merkleTree->verify_merkle_commit(par2);

    std::vector<std::size_t> query_index=merkleTree->find_merkle_path_only_index(2*height-1);
    par3 = merkleTree->create_merklePar_of_mat_by_index(value_for_commit,query_index);
    bool suc4 = merkleTree->verify_merkle_commit(par3);

//    assert(merkleTree->allNodes_.size()==2*width-1)
    if(!(suc3&& suc4 && suc2)){
        std::cout<<"Error\n";
        return 0;
    } else{
        std::cout<<"Suc\n";
    }
    return 0;
}