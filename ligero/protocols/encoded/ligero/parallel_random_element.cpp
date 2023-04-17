#include <vector>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/fields/prime_base/fields_64.hpp>

// 128
//void paralle_random_element(std::vector<libff::alt_bn128_Fr> &result,
//                            std::size_t index,
//                            std::size_t limit){
//    for (std::size_t i = index; i < index+limit; i++) {
//        result[i]=libff::alt_bn128_Fr::random_element();
//    }
//}

//64
void paralle_random_element(std::vector<libff::Fields_64> &result,
                            std::size_t index,
                            std::size_t limit){
    for (std::size_t i = index; i < index+limit; i++) {
        result[i]=libff::Fields_64::random_element();
    }
}
