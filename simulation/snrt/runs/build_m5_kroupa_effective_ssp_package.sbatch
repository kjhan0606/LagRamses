#!/usr/bin/env bash
# Build a source-matched Kroupa v5 package for the approved effective-SSP
# comparison. The source archives stay in the private LRD_JWST asset tree.
#SBATCH --job-name=m5-kroupa-pkg
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=none
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-kroupa-source-package-%j.log

set -euo pipefail

module purge
module load python/3.13.9

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly ASSETS=/gpfs/kjhan/LRD_JWST
readonly SOURCE_DIR=$ASSETS/.parsec-sources.9ZIUlJ
readonly SUKHBOLD_ROOT=$ASSETS/external/g2_candidates
readonly PYTHON=/home/kjhan/miniconda3/bin/python3
readonly WORK=$ROOT/.m5-kroupa-source-package-${SLURM_JOB_ID:?}
readonly EXISTING=$ASSETS/.mixed-lowmass.Rbb2fZ/input-v5

test -d "$SOURCE_DIR"
test -d "$SUKHBOLD_ROOT/sukhbold2016_ccsn"
test -x "$PYTHON"
test -d "$EXISTING"
test ! -e "$WORK"
mkdir -- "$WORK"
exec > >(tee -a "$WORK/build.log") 2>&1

export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

printf 'KROUPA_MIXED_SOURCE_BUILD job=%s python=%s\n' "$SLURM_JOB_ID" "$PYTHON"
"$PYTHON" -c 'import sys, numpy; print("PYTHON",sys.version); print("NUMPY",numpy.__version__)'

"$PYTHON" "$ROOT/simulation/snrt/tools/build_parsec_pair_feedback.py" \
  --source-dir "$SOURCE_DIR" \
  --output "$WORK/upper-feedback" \
  --wind-model fixed \
  --metallicity-grid precision_eleven \
  --wind-km-s 1000 \
  --imf-id 1

"$PYTHON" "$ROOT/simulation/snrt/tools/build_parsec_native_sed.py" \
  --source-dir "$SOURCE_DIR" \
  --feedback "$WORK/upper-feedback" \
  --output "$WORK/upper-sed" \
  --escape-fraction 1.0

"$PYTHON" "$ROOT/simulation/snrt/tools/build_parsec_mixed_lowmass.py" \
  --source-dir "$SOURCE_DIR" \
  --feedback "$WORK/upper-feedback" \
  --upper-sed "$WORK/upper-sed" \
  --output "$WORK/input-v5" \
  --agb-wind-km-s 10 \
  --sn-wind-km-s 1000 \
  --sukhbold-root "$SUKHBOLD_ROOT"

"$PYTHON" -B - "$WORK" "$EXISTING" "$ROOT/simulation/snrt/config/fp2_snia_effective_population_runtime_v1.nml" <<'PY'
import hashlib
import json
from pathlib import Path
import re
import sys

work, previous, snia_contract = map(Path, sys.argv[1:])
package = work / "input-v5"
manifest_path = package / "manifest.json"
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
if manifest.get("model") != "parsec_mixed_lowmass_truncated_v1":
    raise SystemExit("unexpected generated stellar source model")
if manifest.get("imf_id") != 1 or manifest.get("imf_denominator") != [0.08, 600]:
    raise SystemExit("generated package is not normalized for the selected Kroupa IMF")

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

for name, expected in manifest["output_sha256"].items():
    if sha(package / name) != expected:
        raise SystemExit(f"generated package hash mismatch: {name}")

source = (package / "source.nml").read_text(encoding="ascii")
history = (package / "history.nml").read_text(encoding="ascii")
if not re.search(r"\bimf_id=1, population_id=0, binary_fraction=0\b", source):
    raise SystemExit("base source metadata does not bind Kroupa single-star evolution")
if not re.search(r"\binput_imf_id=1\b", history):
    raise SystemExit("history metadata does not bind the selected Kroupa IMF")
node_hash = sha(package / "nodes.dat")
match = re.search(r"source_sha256='([0-9a-f]{64})'", source)
if not match or match.group(1) != node_hash:
    raise SystemExit("SED source identity does not match the generated node history")

# The approved effective-SSP path adds the empirical DTD to single-star
# evolution. Record that population binding as a separate runtime overlay.
overlay = re.sub(r"\bpopulation_id=0, binary_fraction=0\b",
                 "population_id=2, binary_fraction=0", source)
if overlay == source or overlay.count("population_id=2, binary_fraction=0") != 1:
    raise SystemExit("could not create the unique effective-SSP population binding")
overlay_path = work / "source-effective-ssp.nml"
overlay_path.write_text(overlay, encoding="ascii")
if sha(package / "nodes.dat") != node_hash:
    raise SystemExit("node history changed while creating runtime binding")

binding = {
    "schema": "m5-effective-ssp-runtime-binding",
    "schema_version": 1,
    "source_package_manifest_sha256": sha(manifest_path),
    "source_package_model": manifest["model"],
    "source_package_imf_id": manifest["imf_id"],
    "source_package_node_history_sha256": node_hash,
    "source_namelist_sha256": sha(package / "source.nml"),
    "effective_ssp_namelist_sha256": sha(overlay_path),
    "population_id_overlay": {"from": 0, "to": 2},
    "binary_fraction": 0.0,
    "snia_runtime_contract_sha256": sha(snia_contract),
    "interpretation": "single-star PARSEC mixed evolution plus the separately approved empirical effective-SSP SNIa DTD",
    "status": "reference_control_not_production_approval",
}
(work / "effective-ssp-binding.json").write_text(
    json.dumps(binding, indent=2, sort_keys=True) + "\n", encoding="utf-8")

for name in ("yields.dat", "nodes.dat"):
    old_hash, new_hash = sha(previous / name), sha(package / name)
    print(f"SOURCE_REWEIGHT_CHECK {name} previous={old_hash} kroupa={new_hash} identical={old_hash == new_hash}")
print("KROUPA_EFFECTIVE_SSP_PACKAGE_PASS", json.dumps(binding, sort_keys=True))
PY

sha256sum "$WORK/input-v5/manifest.json" "$WORK/input-v5/yields.dat" \
  "$WORK/input-v5/history.nml" "$WORK/input-v5/nodes.dat" \
  "$WORK/input-v5/source.nml" "$WORK/input-v5/enrichment.nml" \
  "$WORK/source-effective-ssp.nml" "$WORK/effective-ssp-binding.json" \
  > "$WORK/package.sha256"
printf 'KROUPA_SOURCE_PACKAGE_BUILD_COMPLETE path=%s\n' "$WORK"
