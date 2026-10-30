/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, and distribute this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 */

#include "node_storage.hpp"

// The storage itself is header-only: it is the implicit-heap layout the builder writes by index and
// the traversals read by index, and it carries no out-of-line behaviour of its own. This translation
// unit stays so the build (and any future explicit-layout code) has somewhere to live.
namespace bvh::nodes {

static_assert(sizeof(Node) == 2 * sizeof(float3), "a node is two bounds vectors");

} // namespace bvh::nodes
