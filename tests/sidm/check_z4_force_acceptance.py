#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def failed_checks(result: dict[str, object]) -> list[str]:
    checks = result["resolution_acceptance_pass"]
    if not isinstance(checks, dict):
        raise RuntimeError("resolution acceptance status is not a mapping")
    failed: list[str] = []
    for scope, values in checks.items():
        if not isinstance(values, dict):
            raise RuntimeError(f"resolution acceptance scope {scope} is not a mapping")
        for name, passed in values.items():
            if passed is not True:
                failed.append(f"{scope}.{name}")
    return failed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result", type=Path)
    args = parser.parse_args()
    result = json.loads(args.result.read_text())
    failed = failed_checks(result)
    if failed:
        raise SystemExit("z~4 force-resolution acceptance failed: " + ", ".join(failed))
    print("Z4_FORCE_RESOLUTION_ACCEPTANCE_RESULT=passed")


if __name__ == "__main__":
    main()
