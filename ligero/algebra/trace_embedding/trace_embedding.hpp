/**@file
 *****************************************************************************
 Trace embedding class definition
 *****************************************************************************
 * @author     This file is part of ligero (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef ligero_ALGEBRA_TRACE_EMBEDDING_TRACE_EMBEDDING_HPP_
#define ligero_ALGEBRA_TRACE_EMBEDDING_TRACE_EMBEDDING_HPP_

#include "ligero/algebra/field_subset/field_subset.hpp"
#include "ligero/algebra/trace_embedding/bivariate_embedding.hpp"
#include "ligero/algebra/trace_embedding/successor_ordering.hpp"

namespace ligero {

template<typename FieldT>
class trace_embedding {
protected:
    bivariate_embedding<FieldT> bivariate_embedding_;
    successor_ordering<FieldT> successor_ordering_;
public:
    trace_embedding<FieldT>() {};
    trace_embedding<FieldT>(field_subset<FieldT> &H, field_subset<FieldT> &row_domain, field_subset<FieldT> &col_domain);
    bivariate_embedding<FieldT> bivariate_embedding() const;
    successor_ordering<FieldT> successor_ordering() const;
};

} // namespace ligero

#include "ligero/algebra/trace_embedding/trace_embedding.tcc"

#endif // ligero_ALGEBRA_TRACE_EMBEDDING_TRACE_EMBEDDING_HPP_
