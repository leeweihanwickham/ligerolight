/** @file
 *****************************************************************************

 Declaration of interfaces for a Merkle tree.

 *****************************************************************************
 * @author     This file is part of libiop, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MERKLE_TREE_HPP_
#define MERKLE_TREE_HPP_

#include <map>
#include <vector>

#include <libff/common/utils.hpp>

namespace ligero {

/**
 * A Merkle tree is maintained as two maps:
 * - a map from addresses to values, and
 * - a map from addresses to hashes.
 *
 * The second map maintains the intermediate hashes of a Merkle tree
 * built atop the values currently stored in the tree (the
 * implementation admits a very efficient support for sparse
 * trees). Besides offering methods to load and store values, the
 * class offers methods to retrieve the root of the Merkle tree and to
 * obtain the authentication paths for (the value at) a given address.
 */

typedef libff::bit_vector merkle_authentication_node;
typedef std::vector<merkle_authentication_node> merkle_authentication_path;

} // ligero


#endif // MERKLE_TREE_HPP_
