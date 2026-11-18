#pragma once

// TaskFunc: fixed-capacity owning callable wrapper.  Requires C++20.
//
//   - No heap allocation: callable lives in an inline byte buffer (SBO).
//   - No virtual functions: type erasure via three plain function pointers.
//   - No C++ standard library headers: traits, forwarding, size_t are hand-rolled;
//     placement new and declval are declared here (freestanding-safe).
//
// C++20 P0593 (implicit object creation) makes placement-new on the byte buffer
// well-defined; no std::launder is needed because every access obtains a fresh
// pointer from the buffer's address, which equals the constructed object's address.
//
// Callables must be nothrow move constructible; enforced at compile time.
// Destructors that throw will terminate, as usual.

// Placement new: provide our own only when <new> has not been pulled in.
// libstdc++ guards <new> with _NEW; libc++ with _LIBCPP_NEW; MSVC with _NEW.
// #if !defined(_NEW) && !defined(_LIBCPP_NEW) && !defined(TF_NO_PLACEMENT_NEW)
// inline void* operator new(decltype(sizeof(0)), void* p) noexcept { return p; }
// #endif

namespace tf {

using TaskFuncSizeType = decltype(sizeof(0));

template <typename A, typename B> struct TaskFuncIsSame { static constexpr bool value = false; };
template <typename A> struct TaskFuncIsSame<A, A> { static constexpr bool value = true; };

template <typename T> struct TaskFuncRemoveRef { using type = T; };
template <typename T> struct TaskFuncRemoveRef<T&> { using type = T; };
template <typename T> struct TaskFuncRemoveRef<T&&> { using type = T; };

template <typename T> struct TaskFuncRemoveCV { using type = T; };
template <typename T> struct TaskFuncRemoveCV<const T> { using type = T; };
template <typename T> struct TaskFuncRemoveCV<volatile T> { using type = T; };
template <typename T> struct TaskFuncRemoveCV<const volatile T> { using type = T; };

// Freestanding declval: no definition, only valid in unevaluated contexts.
template <typename T> T&& taskfunc_declval() noexcept;

[[noreturn]] inline void taskfunc_abort() noexcept {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_trap();
#elif defined(_MSC_VER)
    __debugbreak();
    __assume(0);
#else
    for (;;) {}
#endif
}

template <typename Signature,
          TaskFuncSizeType Size = 32,
          TaskFuncSizeType Alignment = alignof(void*)>
class TaskFunc;

template <typename R, typename... Args,
          TaskFuncSizeType Size, TaskFuncSizeType Alignment>
class TaskFunc<R(Args...), Size, Alignment> {
    // static_assert(!std::is_array_v<R> && !std::is_function_v<R>,
    //               "TaskFunc cannot return an array or function type"); 不好的地方是引入了std 不freestanding
    using InvokeFn  = R (*)(void*, Args&&...);
    using DestroyFn = void (*)(void*);
    using MoveFn    = void (*)(void*, void*); // 擦除后会破坏Destroy自动析构 Move也会被破坏 所以需要如此储存 不需要手动注册

    static_assert(Size >= 1, "TaskFunc storage size must be at least 1");
    static_assert(Alignment >= 1 && (Alignment & (Alignment - 1)) == 0,
                  "TaskFunc alignment must be a power of two");

    InvokeFn  invoke_  = nullptr;
    DestroyFn destroy_ = nullptr;
    MoveFn    move_    = nullptr;
    alignas(Alignment) unsigned char storage_[Size];

    template <typename F>
    using Stored = typename TaskFuncRemoveCV<typename TaskFuncRemoveRef<F>::type>::type;

public:
    TaskFunc() noexcept = default;
    TaskFunc(const TaskFunc&) = delete;
    TaskFunc& operator=(const TaskFunc&) = delete;

    template <typename F>
        requires (!TaskFuncIsSame<Stored<F>, TaskFunc>::value)
              && requires(Stored<F>& s, Args&&... args) {
                     s(static_cast<Args&&>(args)...);
                 }
    TaskFunc(F&& f) {
        using S = Stored<F>;
        static_assert(sizeof(S) <= Size,
                      "callable too large for TaskFunc storage");
        static_assert(alignof(S) <= Alignment,
                      "callable over-aligned for TaskFunc storage");
        static_assert(noexcept(S(taskfunc_declval<S&&>())),
                      "callable must be nothrow move constructible");

        ::new (static_cast<void*>(storage_)) S(static_cast<F&&>(f));
        invoke_ = [](void* p, Args&&... args) -> R {
            return (*static_cast<S*>(p))(static_cast<Args&&>(args)...);
        };
        destroy_ = [](void* p) { static_cast<S*>(p)->~S(); };
        move_ = [](void* dst, void* src) {
            ::new (dst) S(static_cast<S&&>(*static_cast<S*>(src)));
        };
    }

    TaskFunc(TaskFunc&& other) noexcept
        : invoke_(other.invoke_), destroy_(other.destroy_), move_(other.move_) {
        if (move_) {
            move_(storage_, other.storage_);
            other.destroy_(other.storage_);
        }
        other.invoke_  = nullptr;
        other.destroy_ = nullptr;
        other.move_    = nullptr;
    }

    TaskFunc& operator=(TaskFunc&& other) noexcept {
        if (this != &other) {
            reset();
            if (other.move_) {
                invoke_  = other.invoke_;
                destroy_ = other.destroy_;
                move_    = other.move_;
                move_(storage_, other.storage_);
                other.destroy_(other.storage_);
            }
            other.invoke_  = nullptr;
            other.destroy_ = nullptr;
            other.move_    = nullptr;
        }
        return *this;
    }

    ~TaskFunc() { reset(); }

    void reset() noexcept {
        if (destroy_) {
            destroy_(storage_);
        }
        invoke_  = nullptr;
        destroy_ = nullptr;
        move_    = nullptr;
    }

    void swap(TaskFunc& other) noexcept {
        if (this == &other) {
            return;
        }
        TaskFunc tmp(static_cast<TaskFunc&&>(*this));
        *this = static_cast<TaskFunc&&>(other);
        other = static_cast<TaskFunc&&>(tmp);
    }

    explicit operator bool() const noexcept { return invoke_ != nullptr; }

    R operator()(Args... args) const {
        if (!invoke_) {
            taskfunc_abort();
        }
        return invoke_(const_cast<void*>(static_cast<const void*>(storage_)),
                       static_cast<Args&&>(args)...);
    }
};

template <typename Sig, TaskFuncSizeType S, TaskFuncSizeType A>
void swap(TaskFunc<Sig, S, A>& a, TaskFunc<Sig, S, A>& b) noexcept {
    a.swap(b);
}

} // namespace tf
