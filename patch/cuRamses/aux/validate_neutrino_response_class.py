#!/usr/bin/env python3
"""Compare the production CAMB neutrino-response table with CLASS.

The compared quantity is R_nu=T_nu/T_cb, where CLASS ``cb`` is the
density-weighted CDM+baryon transfer.  The gate is restricted to the
production analysis domain z<=2 and 1e-3<=k<=0.5 h/Mpc.  CAMB and CLASS use
independent Einstein--Boltzmann implementations, so agreement establishes a
Level-3 independent linear-response comparison.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np


DEFAULT_TABLE = Path(
    "/gpfs/kjhan/Hydro/DE_nonstd/DMO_production_L512_N1024_20260729/"
    "tables/neutrino_mnu006.dat"
)
DEFAULT_CLASS = Path("/home/kjhan/BACKUP/lagCAMB_validation/class_public/class")


def read_table(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values: list[float] = []
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("#") or not line.strip():
            continue
        values.extend(float(value) for value in line.split())
    nk, na = int(values[0]), int(values[1])
    cursor = 2
    k = np.asarray(values[cursor : cursor + nk])
    cursor += nk
    a = np.asarray(values[cursor : cursor + na])
    cursor += na
    ratio = np.asarray(values[cursor:]).reshape(nk, na)
    if ratio.shape != (nk, na):
        raise ValueError(f"bad table shape {ratio.shape}, expected {(nk, na)}")
    return k, a, ratio


def class_ini(root: Path, redshifts: list[float]) -> str:
    omega_nu = 0.001409
    return f"""h = 0.6766
Omega_b = 0.049
Omega_cdm = {0.3111 - 0.049 - omega_nu:.12f}
N_ur = 2.0328
N_ncdm = 1
m_ncdm = 0.06
T_ncdm = 0.71611
YHe = 0.245
A_s = 2.1e-9
n_s = 0.9665
output = mTk
P_k_max_h/Mpc = 12.0
z_pk = {','.join(str(value) for value in redshifts)}
ncdm_fluid_approximation = 3
l_max_ncdm = 50
tol_ncdm_bg = 1.e-10
tol_ncdm_synchronous = 1.e-10
tol_perturb_integration = 1.e-6
root = {root}
headers = yes
format = class
overwrite_root = yes
"""


def load_class_transfers(directory: Path) -> list[dict[str, np.ndarray]]:
    files = sorted(directory.glob("class_*tk.dat"))
    transfers: list[dict[str, np.ndarray]] = []
    for path in files:
        comment_lines = [
            line for line in path.read_text().splitlines() if line.startswith("#")
        ]
        column_header = next(line for line in comment_lines if "d_cdm" in line)
        names = re.findall(r"\d+:([^\s]+(?:\s+\(h/Mpc\))?)", column_header)
        data = np.loadtxt(path)
        if len(names) != data.shape[1]:
            raise ValueError(
                f"cannot parse CLASS header in {path}: names={names}, comments={comment_lines}"
            )
        transfers.append({name: data[:, index] for index, name in enumerate(names)})
    return transfers


def interpolate_table(
    k_table: np.ndarray,
    a_table: np.ndarray,
    ratio_table: np.ndarray,
    k: np.ndarray,
    a: float,
) -> np.ndarray:
    ia = int(np.searchsorted(a_table, a) - 1)
    ia = max(0, min(ia, len(a_table) - 2))
    ta = (np.log(a) - np.log(a_table[ia])) / (
        np.log(a_table[ia + 1]) - np.log(a_table[ia])
    )
    lower = np.interp(np.log(k), np.log(k_table), ratio_table[:, ia])
    upper = np.interp(np.log(k), np.log(k_table), ratio_table[:, ia + 1])
    return (1.0 - ta) * lower + ta * upper


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--class-executable", type=Path, default=DEFAULT_CLASS)
    parser.add_argument(
        "--json", type=Path, default=Path(__file__).with_name("NEUTRINO_CLASS_LEVEL3_AUDIT.json")
    )
    parser.add_argument(
        "--max-median-relative-error", type=float, default=0.02,
        help="per-redshift median tolerance for R_nu itself",
    )
    parser.add_argument(
        "--max-source-relative-error", type=float, default=5.0e-5,
        help="maximum tolerance for the applied factor 1+(Omega_nu/Omega_cb)R_nu",
    )
    args = parser.parse_args()

    redshifts = [0.0, 0.5, 1.0, 2.0]
    k_table, a_table, ratio_table = read_table(args.table)
    class_root = args.class_executable.resolve().parent
    with tempfile.TemporaryDirectory(prefix="neutrino_class_") as temporary:
        work = Path(temporary)
        ini = work / "neutrino.ini"
        ini.write_text(class_ini(work / "class_", redshifts))
        run = subprocess.run(
            [str(args.class_executable), str(ini)],
            cwd=class_root,
            check=True,
            capture_output=True,
            text=True,
        )
        transfers = load_class_transfers(work)
        if len(transfers) != len(redshifts):
            raise RuntimeError(
                f"CLASS wrote {len(transfers)} transfer files for {len(redshifts)} redshifts"
            )

        # CLASS writes the z=0 file without a z suffix and then follows the
        # z_pk order for the numbered files.  Identify each file through its
        # growth amplitude: |d_cdm| decreases monotonically with redshift.
        transfers.sort(
            key=lambda item: float(np.median(np.abs(item["d_cdm"]))), reverse=True
        )
        records: list[dict[str, float]] = []
        all_errors: list[np.ndarray] = []
        all_source_errors: list[np.ndarray] = []
        omega_b = 0.049
        omega_nu = 0.001408659
        omega_cb = 0.3111 - omega_nu
        omega_cdm = omega_cb - omega_b
        for redshift, transfer in zip(redshifts, transfers, strict=True):
            k = transfer["k"] if "k" in transfer else transfer["k (h/Mpc)"]
            d_cb = (omega_cdm * transfer["d_cdm"] + omega_b * transfer["d_b"]) / (
                omega_cdm + omega_b
            )
            ratio_class = transfer["d_ncdm[0]"] / d_cb
            domain = (k >= 1.0e-3) & (k <= 0.5)
            k_domain = k[domain]
            ratio_camb = interpolate_table(
                k_table, a_table, ratio_table, k_domain, 1.0 / (1.0 + redshift)
            )
            error = np.abs(ratio_class[domain] / ratio_camb - 1.0)
            source_camb = 1.0 + (omega_nu / omega_cb) * ratio_camb
            source_class = 1.0 + (omega_nu / omega_cb) * ratio_class[domain]
            source_error = np.abs(source_class / source_camb - 1.0)
            all_errors.append(error)
            all_source_errors.append(source_error)
            records.append(
                {
                    "redshift": redshift,
                    "points": int(error.size),
                    "median_relative_error": float(np.median(error)),
                    "p95_relative_error": float(np.quantile(error, 0.95)),
                    "max_relative_error": float(np.max(error)),
                    "max_applied_source_relative_error": float(np.max(source_error)),
                }
            )

    maximum = float(max(np.max(error) for error in all_errors))
    maximum_source = float(max(np.max(error) for error in all_source_errors))
    maximum_median = float(max(record["median_relative_error"] for record in records))
    accurate_marker = "accurate_massive_neutrino_transfers=True" in args.table.read_text()
    result = {
        "model": "massive-neutrino linear response and F6+neutrino inherited component",
        "quantity": "R_nu(k,a)=T_nu/T_cb",
        "domain": {"k_hmpc": [1.0e-3, 0.5], "redshifts": redshifts},
        "camb_table": str(args.table.resolve()),
        "camb_table_sha256": hashlib.sha256(args.table.read_bytes()).hexdigest(),
        "class_executable": str(args.class_executable.resolve()),
        "class_executable_sha256": hashlib.sha256(args.class_executable.read_bytes()).hexdigest(),
        "class_stdout_tail": run.stdout.splitlines()[-20:],
        "records": records,
        "accuracy_contract": {
            "maximum_median_R_nu_relative_error": args.max_median_relative_error,
            "maximum_applied_source_relative_error": args.max_source_relative_error,
            "rationale": (
                "R_nu approaches zero at high k, so its pointwise relative error is ill-conditioned; "
                "the full physical source factor is bounded over the complete domain."
            ),
        },
        "accurate_massive_neutrino_transfers_marker": accurate_marker,
        "maximum_median_relative_error": maximum_median,
        "maximum_relative_error": maximum,
        "maximum_applied_source_relative_error": maximum_source,
        "passed": (
            accurate_marker
            and maximum_median < args.max_median_relative_error
            and maximum_source < args.max_source_relative_error
        ),
    }
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.json.write_text(payload)
    print(payload, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
