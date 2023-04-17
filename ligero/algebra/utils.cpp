#include <iostream>
#include "ligero/algebra/utils.hpp"

namespace ligero {

size_t gcd(const size_t a, const size_t b)
{
    return b == 0 ? a : gcd(b, a % b);
}

}


