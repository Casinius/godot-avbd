/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 */

#ifndef AVBD_LIST_RANGE_HPP
#define AVBD_LIST_RANGE_HPP

#include <ranges>

namespace avbd {

// A forward range over an intrusive singly-linked list: Node must expose a public
// `Node *next`. Empty iteration = null head; the end iterator is a null-node iterator,
// so begin()/end() compare by pointer and the range is trivially multi-pass.
template <typename Node>
class next_range {
public:
    class iterator {
    public:
        using value_type = Node *;
        using difference_type = std::ptrdiff_t;
        using iterator_concept = std::forward_iterator_tag;

        iterator() = default;
        explicit iterator(Node *node) : node_(node) {}

        Node *operator*() const { return node_; }
        iterator &operator++() { node_ = node_->next; return *this; }
        iterator operator++(int) { auto tmp = *this; ++*this; return tmp; }
        bool operator==(const iterator &) const = default;

    private:
        Node *node_ = nullptr;
    };

    explicit next_range(Node *head) : head_(head) {}

    iterator begin() const { return iterator{head_}; }
    iterator end() const { return iterator{}; }

private:
    Node *head_;
};

static_assert(std::ranges::forward_range<next_range<struct Rigid>>);

} // namespace avbd

#endif // AVBD_LIST_RANGE_HPP
