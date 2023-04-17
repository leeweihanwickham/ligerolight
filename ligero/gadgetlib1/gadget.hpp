/** @file
 *****************************************************************************
 * @author     This file is part of libiop, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef GADGET_HPP_
#define GADGET_HPP_

#include <ligero/gadgetlib1/protoboard.hpp>

namespace ligero {

template<typename FieldT>
class gadget {
protected:
    protoboard<FieldT> &pb;
    const std::string annotation_prefix;
public:
    gadget(protoboard<FieldT> &pb, const std::string &annotation_prefix="");
};

} // ligero
#include <ligero/gadgetlib1/gadget.tcc>

#endif // GADGET_HPP_

