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

#pragma once
#include <cstddef>

namespace avbd {

// Forward declaration - implemented by classes that want to be observed
template<typename T>
class ObserverList;

template<typename T>
class ObservePtr {
public:
    // Default constructor
    ObservePtr() = default;

    // Construct from raw pointer - registers as observer automatically
    ObservePtr(T* p) : m_ptr(p) {
        if (p) {
            p->addObserver(this);
        }
    }

    // Copy constructor
    ObservePtr(const ObservePtr& other) : m_ptr(other.m_ptr) {
        if (m_ptr) {
            m_ptr->addObserver(this);
        }
    }

    // Copy assignment
    ObservePtr& operator=(const ObservePtr& other) {
        if (this != &other) {
            if (m_ptr) {
                m_ptr->removeObserver(this);
            }
            m_ptr = other.m_ptr;
            if (m_ptr) {
                m_ptr->addObserver(this);
            }
        }
        return *this;
    }

    // Move constructor: replace the source node in the observed object's list.
    ObservePtr(ObservePtr&& other) noexcept
        : m_ptr(other.m_ptr), m_prev(other.m_prev), m_next(other.m_next) {
        if (m_ptr)
            m_ptr->replaceObserver(&other, this);
        other.m_ptr = nullptr;
        other.m_prev = nullptr;
        other.m_next = nullptr;
    }

    // Move assignment
    ObservePtr& operator=(ObservePtr&& other) noexcept {
        if (this != &other) {
            if (m_ptr)
                m_ptr->removeObserver(this);
            m_ptr = other.m_ptr;
            m_prev = other.m_prev;
            m_next = other.m_next;
            if (m_ptr)
                m_ptr->replaceObserver(&other, this);
            other.m_ptr = nullptr;
            other.m_prev = nullptr;
            other.m_next = nullptr;
        }
        return *this;
    }

    // Destructor - deregisters automatically
    ~ObservePtr() {
        if (m_ptr) {
            m_ptr->removeObserver(this);
        }
    }

    // Get the raw pointer
    [[nodiscard]] T* get() const { return m_ptr; }

    // Pointer operators
    [[nodiscard]] T* operator->() const { return m_ptr; }
    [[nodiscard]] T& operator*() const { return *m_ptr; }
    
    // Conversion to bool
    [[nodiscard]] explicit operator bool() const { return m_ptr != nullptr; }

    // Called by T while walking its observer list.
    void invalidate() {
        m_ptr = nullptr;
        m_prev = nullptr;
        m_next = nullptr;
    }

    [[nodiscard]] bool is_valid() const { return m_ptr != nullptr; }

private:
    template<typename U>
    friend class ObservePtr;
    friend T;

    // Linked list pointers for the observer chain
    T* m_ptr = nullptr;
    ObservePtr<T>* m_prev = nullptr;
    ObservePtr<T>* m_next = nullptr;
};

} // namespace avbd
