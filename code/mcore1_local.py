"""
Local MCORE-1 metrical-tree builder and checker.

Used when ``vendor/mcore-1`` is not installed. API mirrors the intended upstream
``mcore1.check_tree`` contract: post-order validation for CONSERVATION, OVERFLOW,
and EMPTY_CONSTITUENT.
"""

from __future__ import annotations

from mcore1_types import TreeNode


def build_metrical_tree(trits: list[int]) -> TreeNode:
    """Group emitted trits bottom-up in triples (ternary constituency tree)."""
    leaves = [TreeNode(weight=int(t), children=[]) for t in trits]
    if not leaves:
        return TreeNode(weight=0, children=[], empty=True)

    pad = (3 - len(leaves) % 3) % 3
    for _ in range(pad):
        leaves.append(TreeNode(weight=0, children=[], empty=True))

    level = leaves
    while len(level) > 1:
        nxt: list[TreeNode] = []
        for i in range(0, len(level), 3):
            group = level[i : i + 3]
            s = sum(n.weight for n in group) + 3 * sum(n.overflow for n in group)
            nxt.append(TreeNode(weight=s % 3, overflow=s // 3, children=group))
        level = nxt
    return level[0]


def check_tree(node: TreeNode) -> list[str]:
    """Post-order validation; returns human-readable violation codes."""
    errors: list[str] = []

    if node.empty:
        if node.children:
            errors.append("EMPTY_CONSTITUENT")
        return errors

    if not node.children:
        if node.weight not in (0, 1, 2):
            errors.append(f"INVALID_LEAF_WEIGHT:{node.weight}")
        return errors

    for child in node.children:
        errors.extend(check_tree(child))

    child_sum = sum(c.weight for c in node.children)
    if node.weight != child_sum % 3:
        errors.append(
            f"CONSERVATION: parent={node.weight} children_sum%3={child_sum % 3}"
        )

    expected_overflow = sum(c.overflow for c in node.children) + child_sum // 3
    if node.overflow != expected_overflow:
        errors.append(
            f"OVERFLOW: parent={node.overflow} expected={expected_overflow}"
        )

    return errors
