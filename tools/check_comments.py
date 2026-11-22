#!/usr/bin/env python3
"""Fail if any explanatory comment in src/ is longer than two lines.

The leading copyright/license block is exempt, as is the third-party linalg.h.
Run this after any comment cleanup: python3 tools/check_comments.py
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
EXEMPT = {"linalg.h"}
LICENSE_LINES = 14


def is_comment(line: str) -> bool:
    s = line.lstrip()
    return s.startswith("//") or s.startswith("/*") or s.startswith("*") or s.startswith("*/")


violations = []
for path in sorted((ROOT / "src").rglob("*")):
    if path.suffix not in {".cpp", ".h", ".hpp"} or path.name in EXEMPT:
        continue
    lines = path.read_text(errors="replace").splitlines()
    i = LICENSE_LINES
    while i < len(lines):
        if is_comment(lines[i]):
            j = i
            while j < len(lines) and is_comment(lines[j]):
                j += 1
            if j - i > 2:
                violations.append(f"{path.relative_to(ROOT)}:{i + 1} ({j - i} lines)")
            i = j
        else:
            i += 1

if violations:
    print("comments longer than two lines:")
    print("\n".join("  " + v for v in violations))
    sys.exit(1)
print(f"OK: all comments in src/ are at most two lines ({len(violations)} violations)")
