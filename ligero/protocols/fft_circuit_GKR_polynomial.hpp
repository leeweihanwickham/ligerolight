//
// Created by wx on 22-11-26.
//

#ifndef LIGERO_FFT_CIRCUIT_GKR_POLYNOMIAL_HPP
#define LIGERO_FFT_CIRCUIT_GKR_POLYNOMIAL_HPP

/** @brief several univariate polynomial definitions.
    */
#include <vector>
namespace ligero{

template<typename FieldT>
class linear_poly;
/**Construct a univariate cubic polynomial of f(x) = ax^3 + bx^2 + cx + d
        */
template<typename FieldT>
class cubic_poly
{
public:
    FieldT a, b, c, d;
    cubic_poly();
    cubic_poly(const FieldT&, const FieldT&, const FieldT&, const FieldT&);
    cubic_poly operator + (const cubic_poly<FieldT> &) const;
    FieldT eval(const FieldT &) const;
};

/**Construct a univariate quadratic polynomial of f(x) = ax^2 + bx + c
        */
template<typename FieldT>
class quadratic_poly
{
public:
    FieldT a, b, c;
    quadratic_poly();
    quadratic_poly(const FieldT&, const FieldT&, const FieldT&);
    quadratic_poly<FieldT> operator + (const quadratic_poly<FieldT> &) const;
    cubic_poly<FieldT> operator * (const linear_poly<FieldT> &) const;
    FieldT eval(const FieldT &) const;
};

/**Construct a univariate linear polynomial of f(x) = ax + b
        */
template<typename FieldT>
class linear_poly
{
public:
    FieldT a, b;
    linear_poly();
    linear_poly(const FieldT &, const FieldT &);
    linear_poly(const FieldT &);
    linear_poly<FieldT> operator + (const linear_poly<FieldT> &) const;
    quadratic_poly<FieldT> operator * (const linear_poly<FieldT> &) const;
    FieldT eval(const FieldT &) const;
};



/**Construct a univariate quintuple polynomial of f(x) = ax^4 + bx^3 + cx^2 + dx + e
        */
template<typename FieldT>
class quadruple_poly
{
public:
    FieldT a, b, c, d, e;
    quadruple_poly();
    quadruple_poly(const FieldT&, const FieldT&, const FieldT&, const FieldT&, const FieldT&);
    quadruple_poly<FieldT> operator + (const quadruple_poly<FieldT> &) const;
    FieldT eval(const FieldT &) const;
};

/**Construct a univariate quintuple polynomial of f(x) = ax^5 + bx^4 + cx^3 + dx^2 + ex + f
        */
template<typename FieldT>
class quintuple_poly
{
public:
    FieldT a, b, c, d, e, f;
    quintuple_poly();
    quintuple_poly(const FieldT&, const FieldT&, const FieldT&, const FieldT&, const FieldT&, const FieldT&);
    quintuple_poly<FieldT> operator + (const quintuple_poly<FieldT> &) const;
    FieldT eval(const FieldT &) const;
};

template<typename FieldT>
quintuple_poly<FieldT>::quintuple_poly(){}

template<typename FieldT>
quintuple_poly<FieldT>::quintuple_poly(const FieldT& aa, const FieldT& bb, const FieldT& cc, const FieldT& dd, const FieldT& ee, const FieldT& ff)
{
    a = aa;
    b = bb;
    c = cc;
    d = dd;
    e = ee;
    f = ff;
}

template<typename FieldT>
quintuple_poly<FieldT> quintuple_poly<FieldT>::operator + (const quintuple_poly &x) const
{
    return quintuple_poly<FieldT>(a + x.a, b + x.b, c + x.c, d + x.d, e + x.e, f + x.f);
}

template<typename FieldT>
FieldT quintuple_poly<FieldT>::eval(const FieldT &x) const
{
    return (((((a * x) + b) * x + c) * x + d) * x + e) * x + f;
}

/* quadruple_poly implementation */
template<typename FieldT>
quadruple_poly<FieldT>::quadruple_poly(){}

template<typename FieldT>
quadruple_poly<FieldT>::quadruple_poly(const FieldT& aa, const FieldT& bb, const FieldT& cc, const FieldT& dd, const FieldT& ee)
{
    a = aa;
    b = bb;
    c = cc;
    d = dd;
    e = ee;
}

template<typename FieldT>
quadruple_poly<FieldT> quadruple_poly<FieldT>::operator + (const quadruple_poly<FieldT> &x) const
{
    return quadruple_poly(a + x.a, b + x.b, c + x.c, d + x.d, e + x.e);
}

template<typename FieldT>
FieldT quadruple_poly<FieldT>::eval(const FieldT &x) const
{
    return ((((a * x) + b) * x + c) * x + d) * x + e;
}

/* cubic_poly implementation */
template<typename FieldT>
cubic_poly<FieldT>::cubic_poly(){}
template<typename FieldT>
cubic_poly<FieldT>::cubic_poly(const FieldT& aa, const FieldT& bb, const FieldT& cc, const FieldT& dd)
{
    a = aa;
    b = bb;
    c = cc;
    d = dd;
}

template<typename FieldT>
cubic_poly<FieldT> cubic_poly<FieldT>::operator + (const cubic_poly<FieldT> &x) const
{
    return cubic_poly(a + x.a, b + x.b, c + x.c, d + x.d);
}

template<typename FieldT>
FieldT cubic_poly<FieldT>::eval(const FieldT &x) const
{
    return (((a * x) + b) * x + c) * x + d;
}

/* quadratic_poly implementation */
template<typename FieldT>
quadratic_poly<FieldT>::quadratic_poly(){}
template<typename FieldT>
quadratic_poly<FieldT>::quadratic_poly(const FieldT& aa, const FieldT& bb, const FieldT& cc)
{
    a = aa;
    b = bb;
    c = cc;
}

template<typename FieldT>
quadratic_poly<FieldT> quadratic_poly<FieldT>::operator + (const quadratic_poly<FieldT> &x) const
{
    return quadratic_poly<FieldT>(a + x.a, b + x.b, c + x.c);
}

template<typename FieldT>
cubic_poly<FieldT> quadratic_poly<FieldT>::operator * (const linear_poly<FieldT> &x) const
{
    return cubic_poly<FieldT>(a * x.a, a * x.b + b * x.a, b * x.b + c * x.a, c * x.b);
}

template<typename FieldT>
FieldT quadratic_poly<FieldT>::eval(const FieldT &x) const
{
    return ((a * x) + b) * x + c;
}




/* linear_poly implementation */
template<typename FieldT>
linear_poly<FieldT>::linear_poly(){}

template<typename FieldT>
linear_poly<FieldT>::linear_poly(const FieldT& aa, const FieldT& bb)
{
    a = aa;
    b = bb;
}

template<typename FieldT>
linear_poly<FieldT>::linear_poly(const FieldT &x)
{
    a = FieldT(0);
    b = x;
}

template<typename FieldT>
linear_poly<FieldT> linear_poly<FieldT>::operator + (const linear_poly<FieldT> &x) const
{
    return linear_poly<FieldT>(a + x.a, b + x.b);
}

template<typename FieldT>
quadratic_poly<FieldT> linear_poly<FieldT>::operator * (const linear_poly<FieldT> &x) const
{
    return quadratic_poly<FieldT>(a * x.a, a * x.b + b * x.a, b * x.b);
}

template<typename FieldT>
FieldT linear_poly<FieldT>::eval(const FieldT &x) const
{
    return a * x + b;
}

}

#endif //LIGERO_FFT_CIRCUIT_GKR_POLYNOMIAL_HPP
