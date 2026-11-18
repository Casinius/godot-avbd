#pragma once
#include "func/taskfunc.hpp"
// ============================================================================
//  task.hpp —— 与 State 完全无关的执行器（任务层的唯一注入口）
// ----------------------------------------------------------------------------
//  【为什么单独一个文件】
//  求解器核心（state.hpp / base_concept.hpp）里**不允许**出现任何线程原语：
//  没有 <thread>、没有 std::atomic、没有 mutex。并发只能通过这里定义的 Executor
//  以参数形式注入，这样：
//    - 求解器的语义与线程数解耦，"串行执行器"就是可证的参考实现；
//    - 换 TBB / 自研调度器 / GPU 后端时，只需要再实现一个 Executor；
//    - 上层想怎么排任务（帧流水线、跨帧重叠）都不必改动求解器一行。
//  tests/run.sh 里有 grep 门禁在执行这条约束。
//
//  【契约（重要，别绕）】
//   1) parallel_for(n, grain, fn) 保证对 [0, n) 的每个 index **恰好一次**。
//      分发是**静态块划分**（块号 = index / grain），调用顺序/线程归属不保证。
//   2) **不可重入**：在 fn 内部再调用 parallel_for 会死锁。求解阶段必须
//      在顶层顺序调用（这也是为什么 solve() 的阶段划分是线性的）。
//   3) fn 拿到的区间是 [begin, end)，永不为空。
//   4) 确定性靠调用方：本文件不保证跨线程数的浮点求和顺序，只保证"每个 index
//      被执行一次"。要位一致就得按 index 分区 / 用 slot（求解器就是这么做的）。
//   5) 异常是关闭的（-fno-exceptions），这里也就不做任何异常传播。
//
//  【性能注记 —— 为什么是 spin-then-park】
//  物理一步里有 iters×4 次 fork-join（chain/32 是 800 次）。用"纯条件变量"唤醒时，
//  一次 fork-join 的固定成本实测 17~20 us（每个 worker 要穿两次 mutex，再加主线程
//  的 done 条件变量），800 次 ≈ 14 ms —— 直接把 4 ms 的活干成 18 ms。
//  这里改成：worker 先自旋若干次（阶段间只隔几十微秒，几乎总能自旋命中），
//  超时才挂起；完成计数用原子，主线程也是先自旋后挂起。
//  典型 fork-join 固定成本据此降到 ~1-2 us 量级（见 `avbd-mo --bench`）。
// ============================================================================

#include <algorithm>
#include <atomic>
#include <concepts>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace task {

// 自旋提示：x86 的 pause / arm 的 yield；其它平台退回让出时间片。
inline void spin_hint() noexcept {
#if defined(__x86_64__) || defined(__i386__)
  __builtin_ia32_pause();
#elif defined(__aarch64__) || defined(__arm__)
  asm volatile("yield" : : : "memory");
#else
  std::this_thread::yield();
#endif
}

// ==================== 0. Executor 概念 ====================
// f(begin, end) 消费一个连续区间；这是所有 phase 的通用形态
// （按 body 分块、按 constraint 分块、按 entry 分块都一样）。
template <typename E>
concept Executor = requires(E &e, std::size_t n, std::size_t grain,
                            tf::TaskFunc<void(std::size_t, std::size_t)> &f) {
  { e.parallel_for(n, grain, f) } noexcept;
  { E::is_serial() } -> std::convertible_to<bool>;
  { e.thread_count() } -> std::convertible_to<unsigned>;
};

// ==================== 1. SerialExecutor：参考实现 ====================
// 不是"退化版"——它就是基准语义。所有并发化后的结果都必须与它逐位相同。
class SerialExecutor {
public:
  static constexpr bool is_serial() noexcept { return true; }
  unsigned              thread_count() const noexcept { return 1; }

  template <typename F> void parallel_for(std::size_t n, std::size_t /*grain*/, F &&fn) noexcept {
    if (n != 0)
      fn(std::size_t(0), n);
  }
};

// ==================== 2. ThreadPoolExecutor：fork-join，spin-then-park ====================
// 为什么是 fork-join 而不是通用任务图：
//   一次物理迭代只有 scatter / gather / primal / dual 四个阶段，阶段之间必须
//   全序（barrier）。任务图能提供的重叠在这里不存在，只会额外增加调度开销。
//
// 为什么用静态块划分 + 原子取块，而不是 dynamic work stealing：
//   stealing 的迁移会让"浮点求和顺序"变成时序的函数；块划分至少让块的边界固定，
//   配合调用方的 per-index 隔离就能做位一致。竞争只在每块取号一次（relaxed）。
class ThreadPoolExecutor {
public:
  explicit ThreadPoolExecutor(unsigned threads = 0) {
    if (threads == 0)
      threads = std::thread::hardware_concurrency();
    if (threads == 0)
      threads = 1;
    workers_.reserve(threads > 0 ? threads - 1 : 0);
    for (unsigned i = 1; i < threads; ++i)
      workers_.emplace_back([this] { worker_loop(); });
  }

  ~ThreadPoolExecutor() {
    stop_.store(true, std::memory_order_release);
    start_cv_.notify_all(); // parked worker 的谓词里有 stop_
    for (auto &t : workers_)
      t.join();
  }

  ThreadPoolExecutor(const ThreadPoolExecutor &)            = delete;
  ThreadPoolExecutor &operator=(const ThreadPoolExecutor &) = delete;

  static constexpr bool is_serial() noexcept { return false; }
  unsigned thread_count() const noexcept { return static_cast<unsigned>(workers_.size()) + 1u; }

  // 见文件头的契约 1)~5)。
  template <typename F> void parallel_for(std::size_t n, std::size_t grain, F &&fn) noexcept {
    if (n == 0)
      return;
    if (grain == 0)
      grain = 1;
    if (workers_.empty() || n <= grain) { // 太小，fork 的开销比工作本身还大
      fn(std::size_t(0), n);
      return;
    }

    const std::size_t chunks = (n + grain - 1) / grain;
    // 只做一次类型擦除（捕获全是引用/标量，且走 tf::TaskFunc 的 SBO，不存在任何堆分配）。
    tf::TaskFunc<void(std::size_t)> job = [&fn, grain, n](std::size_t ci) {
      const std::size_t b = ci * grain;
      fn(b, std::min(b + grain, n));
    };

    publish(&job, chunks);
    if (parked_.load(std::memory_order_acquire) != 0)
      start_cv_.notify_all();

    drain(chunks); // 调用线程也干活，别让 worker 空转等

    // 先自旋等完成（阶段间很密，通常几步就命中），超时才挂起。
    const std::size_t want = workers_.size();
    for (int spin = 0; spin < kSpinIters; ++spin) {
      if (workers_done_.load(std::memory_order_acquire) == want)
        return;
      spin_hint();
    }
    std::unique_lock<std::mutex> lk(m_);
    done_cv_.wait(lk,
                  [this, want] { return workers_done_.load(std::memory_order_acquire) == want; });
  }

private:
  static constexpr int kSpinIters = 4000;

  void publish(tf::TaskFunc<void(std::size_t)> *job, std::size_t chunks) noexcept {
    std::lock_guard<std::mutex> lk(m_);
    job_    = job;
    chunks_ = chunks;
    next_.store(0, std::memory_order_relaxed);
    workers_done_.store(0, std::memory_order_relaxed);
    generation_.fetch_add(1, std::memory_order_release); // release：发布 job_/chunks_
  }

  void drain(std::size_t chunks) noexcept {
    for (;;) {
      const std::size_t ci = next_.fetch_add(1, std::memory_order_relaxed);
      if (ci >= chunks)
        return;
      (*job_)(ci);
    }
  }

  void park(std::size_t seen) noexcept {
    std::unique_lock<std::mutex> lk(m_);
    parked_.fetch_add(1, std::memory_order_seq_cst);
    start_cv_.wait(lk, [this, seen] {
      return stop_.load(std::memory_order_acquire) ||
             generation_.load(std::memory_order_acquire) != seen;
    });
    parked_.fetch_sub(1, std::memory_order_seq_cst);
  }

  void worker_loop() noexcept {
    std::size_t seen = 0;
    for (;;) {
      // —— 等新 generation（或 stop）——
      for (int spin = 0;; ++spin) {
        if (stop_.load(std::memory_order_acquire))
          return;
        const std::size_t g = generation_.load(std::memory_order_acquire);
        if (g != seen) {
          seen = g;
          break;
        }
        if (spin >= kSpinIters) {
          park(seen); // 返回后重新检查 stop_ / generation_
        } else {
          spin_hint();
        }
      }
      if (stop_.load(std::memory_order_acquire))
        return;

      drain(chunks_);
      if (workers_done_.fetch_add(1, std::memory_order_acq_rel) + 1 == workers_.size()) {
        // 只有最后一个完成的 worker 负责叫醒主线程。持锁通知，避免丢唤醒。
        std::lock_guard<std::mutex> lk(m_);
        done_cv_.notify_all();
      }
    }
  }

  std::vector<std::thread> workers_;
  std::mutex               m_;
  std::condition_variable  start_cv_; // 主线程 → worker（parked 时才用）
  std::condition_variable  done_cv_;  // worker → 主线程（自旋超时才用）

  tf::TaskFunc<void(std::size_t)> *job_    = nullptr;
  std::size_t                      chunks_ = 0;

  std::atomic<std::size_t> generation_{0};
  std::atomic<std::size_t> next_{0};
  std::atomic<std::size_t> workers_done_{0};
  std::atomic<std::size_t> parked_{0};
  std::atomic<bool>        stop_{false};
};

} // namespace task
