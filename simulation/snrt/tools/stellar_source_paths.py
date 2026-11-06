"""Resolve explicitly scoped stellar inputs without borrowing another project.

Legacy callers retain their existing defaults. An alternate contract must
declare ``source_root`` or receive it explicitly; relative paths then refer
to that root. Paths are resolved through symlinks before any source is read.
This module changes input routing and provenance, never source physics.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

from adapt_g2_candidate_sources import SourceAdapterError


class SourcePathError(SourceAdapterError):
    """An explicitly scoped source refers outside its declared root."""


def source_path(path, *, root=None):
    path = Path(path)
    if root is not None:
        root = Path(root).resolve()
        if not path.is_absolute():
            path = root / path
        path = path.resolve()
        if not path.is_relative_to(root):
            raise SourcePathError(f"source path escapes declared root {root}: {path}")
        return path
    return path.resolve()


def source_child(base, name, *, root=None):
    """Disallow traversal; symlink targets must stay in the declared source root.

    A scoped view may link to another raw-source directory inside that root.
    Without an explicit root, retain legacy reads through external symlinks.
    """
    base = source_path(base, root=root)
    name = Path(name)
    if name.is_absolute() or ".." in name.parts:
        raise SourcePathError(f"unsafe source manifest filename: {name}")
    return source_path(base / name, root=root)


def contract_root(path, contract, *, root=None, default_path=None):
    path = Path(path).resolve()
    declared = contract.get("source_root")
    if declared is not None:
        declared = Path(declared)
        if not declared.is_absolute():
            declared = path.parent / declared
        declared = declared.resolve()
        if root is not None and Path(root).resolve() != declared:
            raise SourcePathError("supplied source root disagrees with contract source_root")
        root = declared
    if root is None:
        if default_path is None or path != Path(default_path).resolve():
            raise SourcePathError("an alternate source contract requires an explicit source_root")
        return None
    root = Path(root).resolve()
    source_path(path, root=root)
    return root


def file_identity(path, *, root=None):
    path = source_path(path, root=root)
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    return {"path": str(path), "bytes": size, "sha256": digest.hexdigest()}
