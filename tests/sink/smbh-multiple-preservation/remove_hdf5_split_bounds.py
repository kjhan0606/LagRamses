#!/usr/bin/env python3
"""Make a copied checkpoint emulate the legacy, lossy Hilbert-bound format."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    args = parser.parse_args()
    with h5py.File(args.checkpoint, "r+") as stream:
        domain = stream["domain"]
        if int(domain.attrs["bound_key_split_format"][0]) != 1:
            raise ValueError("source checkpoint has no split Hilbert bounds")
        del domain.attrs["bound_key_split_format"]
        del domain["bound_key_hi"]
        del domain["bound_key_lo"]


if __name__ == "__main__":
    main()
