# Taste
- Writes in Chinese and expects replies in Chinese (technical terms/identifiers stay in English). Confidence: 0.6
- Wants performance work to start from real measurements: profile with timers and compare against a reference implementation (e.g. Jolt) under identical conditions (same scene, same physics tick) before proposing optimizations. Confidence: 0.6
- Prefers a written, prioritized optimization plan before code changes: ranked items backed by profiler data, expected gains, and verification gates, rather than piling on more heuristics. Confidence: 0.55
- Drives multi-step work from an explicit todo list: wants the task list created up front and then executed item by item, rather than ad-hoc fixing. Confidence: 0.6
- Treats simulation correctness/stability as outranking performance: a speed win (tighter AABBs, skipped iterations, early exit) that lets dynamic bodies tunnel through static ones is unacceptable, so benchmark gains must be checked against actual scene behaviour, not just ms/step. Confidence: 0.6
- Reports bugs as short, colloquial one-liners in Chinese with no repro steps or stack details, and expects the agent to reproduce them from the shipped demo/test scenes and find the root cause itself. Confidence: 0.55
