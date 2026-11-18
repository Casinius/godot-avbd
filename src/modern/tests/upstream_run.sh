#!/usr/bin/env bash
# ============================================================================
#  tests/run.sh —— 唯一入口：构建 + 结构自检 + 逐位基线比对 + 分层门禁
# ----------------------------------------------------------------------------
#  用法：
#    tests/run.sh            全跑（构建 → 自检 → 基线 → 门禁）
#    tests/run.sh --rebase   确认行为变更后，重新录制基线
#    tests/run.sh --no-build 跳过构建（用现有二进制）
#
#  为什么"基线比对"和"结构自检"要一起跑：
#    基线只能回答"和上次一样吗"，它能被"稳定地错"骗过。自检断言 invariants
#    （dq 一致性、静态 body 冻结、buffers 清空、C==0 时真不动）。两者互补：
#    自检防"错得有理"，基线防"对得更动"。
#
#  分层门禁（step 2 起生效，现在就该是绿的）：
#    求解器核心 state.hpp / base_concept.hpp 不允许出现线程原语。并发只能通过
#    注入的 Executor 进来。这条用 grep 守着，防止"顺手在 State 里加个 mutex"。
# ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."

BIN=build/linux/x86_64/release/avbd-mo
BASELINE=tests/baseline.txt
COLLISION_BASELINE=tests/collision_observable.txt

# --no-build 且 xmake 被配成 debug 时，二进制在另一个目录。
# 不做这个检查的话，脚本会直接报 "missing binary"，掩盖真正原因（模式不对）。
if [[ ! -x $BIN && -x build/linux/x86_64/debug/avbd-mo ]]; then
  BIN=build/linux/x86_64/debug/avbd-mo
fi

REBASE=0
BUILD=1
for arg in "$@"; do
  case "$arg" in
    --rebase)   REBASE=1 ;;
    --no-build) BUILD=0 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done

if [[ $BUILD -eq 1 ]]; then
  echo "== build =="
  # xmake 只按 mtime 判新旧，脚本/编辑器回滚文件内容时不会触发重编。
  # 先 touch 一遍所有源文件，保证二进制一定对应当前源码。
  # 注意：必须用 find 递归，`src/*.hpp` 只覆盖顶层，碰不到 src/phys_core/ 与 src/func/。
  find src -type f \( -name '*.hpp' -o -name '*.cpp' \) -exec touch {} +
  xmake -r
fi

if [[ ! -x $BIN ]]; then
  echo "missing binary: $BIN" >&2
  exit 2
fi

echo
echo "== structural invariants =="
"$BIN" --selftest >/dev/null

echo
echo "== layering gate =="
# 线程原语只允许出现在 task.hpp（执行器实现）里。
# scene/digest 也算核心：它们借用 State 的接口读状态，一旦在那里开线程，
# "State 不知道线程存在"这个性质就没了。
#
# 【别改弱这个门禁】它是静态 grep，不是编译器。曾经因为重构把 state.hpp /
# base_concept.hpp 移到了 src/phys_core/，而这里还写着旧路径，配合
# `2>/dev/null || true` 把“文件不存在”吞掉，结果这条门禁**永远绿**。
# 所以现在：文件缺失 = FAIL，grep 的任何非 0/1 退出码 = FAIL。
GATE_FILES=(src/phys_core/state.hpp src/phys_core/base_concept.hpp src/scene.hpp src/digest.hpp \
             src/collision/shape.hpp src/collision/aabb.hpp src/collision/concepts.hpp src/collision/broadphase.hpp \
             src/collision/narrowphase.hpp src/collision/manifold.hpp src/collision/contact_constraints.hpp \
             src/collision/collision_world.hpp src/collision/collision_digest.hpp \
             src/collision/shape_instance.hpp src/collision/shape_table.hpp \
             src/collision/convex_mesh.hpp src/collision/convex_hull.hpp \
             src/collision/mass_properties.hpp)
for f in "${GATE_FILES[@]}"; do
  if [[ ! -f $f ]]; then
    echo "FAIL: layering gate cannot run, missing $f" >&2
    echo "      (门禁文件被移动/改名了：请同步更新 tests/run.sh 的 GATE_FILES)" >&2
    exit 1
  fi
done
set +e
violations=$(grep -nE '#include <(thread|atomic|mutex|future|condition_variable|semaphore|execution)>|\bstd::(thread|atomic|mutex|jthread|async)\b' \
             "${GATE_FILES[@]}")
grep_rc=$?
set -e
if [[ $grep_rc -gt 1 ]]; then
  echo "FAIL: layering gate grep failed (rc=$grep_rc)" >&2
  exit 1
fi
if [[ -n "$violations" ]]; then
  echo "FAIL: thread primitives leaked out of task.hpp:" >&2
  echo "$violations" >&2
  exit 1
fi
echo "ok: solver core is thread-primitive free (only task.hpp may use them)"

# 碰撞层还多一条：不允许拥有型裸指针 / 哈希容器。碰撞层的 body 引用必须是索引，
# manifold 存储必须是排序 vector；这些规则一旦破掉，确定性契约就守不住。
COLLISION_HEADERS=(src/collision/*.hpp src/collision/bvh/*.hpp)
set +e
collision_raw=$(grep -nE '\bunordered_(map|set)\b|\bRigid\s*\*|\bvoid\s*\*' "${COLLISION_HEADERS[@]}")
grep_rc=$?
set -e
if [[ $grep_rc -gt 1 ]]; then
  echo "FAIL: collision gate grep failed (rc=$grep_rc)" >&2
  exit 1
fi
# 只检查代码行，注释里解释“为什么不用 Rigid*/unordered_map”不算违规。
collision_violations=$(printf '%s\n' "$collision_raw" | grep -vE '^[^:]+:[0-9]+:[[:space:]]*//' || true)
if [[ -n "$collision_violations" ]]; then
  echo "FAIL: collision headers use owning raw pointers / hash containers:" >&2
  echo "$collision_violations" >&2
  exit 1
fi
echo "ok: collision core uses indices + sorted vectors (no owning raw pointers / hash tables)"

echo
if [[ -f $COLLISION_BASELINE ]]; then
  echo "== collision observable digest (warning-only) =="
  if "$BIN" --check-collision "$COLLISION_BASELINE" >/dev/null 2>&1; then
    echo "ok: collision observable digest matches $COLLISION_BASELINE"
  else
    echo "WARN: collision observable digest changed." >&2
    echo "      This is expected when the narrowphase algorithm changes." >&2
    echo "      Review the diff, then run: $BIN --save-collision $COLLISION_BASELINE" >&2
    "$BIN" --check-collision "$COLLISION_BASELINE" 2>&1 | sed 's/^/      /' >&2 || true
  fi
fi

echo
if [[ $REBASE -eq 1 ]]; then
  echo "== rebase baseline =="
  "$BIN" --save "$BASELINE" >/dev/null
  echo "wrote $BASELINE ($(wc -l < "$BASELINE") cases)"
else
  echo "== bit-exact baseline =="
  "$BIN" --check "$BASELINE" >/dev/null

  # 并发正确性：同一张用例表在不同线程数下必须给出逐位相同的摘要。
  # 这是 step 3 起就能成立的性质（primal 逐 body / dual 逐约束都是独立写），
  # 以后任何"顺手并行化"只要引入真实写冲突，这里就会红。
  echo
  echo "== bit-exact across thread counts =="
  for t in 2 3 4 8; do
    "$BIN" --threads "$t" --check "$BASELINE" >/dev/null 2>&1 \
      && echo "ok: threads=$t" \
      || { echo "FAIL: threads=$t diverged from baseline" >&2; exit 1; }
  done

  # 确定性契约（state.hpp 不变量 10）：分区粒度与串行阈值也不许改变结果。
  # grain=1 / min_items=1 是"最碎、最并行"的极端配置。
  echo
  echo "== bit-exact across grain / thresholds =="
  for g in 1 5 16 64; do
    "$BIN" --threads 8 --grain-body "$g" --grain-constraint "$g" --min-items-body 1 \
           --min-items-constraint 1 --check "$BASELINE" >/dev/null 2>&1 \
      && echo "ok: grain=$g min_items=1" \
      || { echo "FAIL: grain=$g diverged from baseline" >&2; exit 1; }
  done
fi

echo
echo "ALL GREEN"
