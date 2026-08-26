#!/bin/bash
#SBATCH --job-name=sidm_zoom_ic
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=2-00:00:00
#SBATCH --output=prepare_zoom_%j.log
#SBATCH --error=prepare_zoom_%j.err

set -euo pipefail

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
parent_run=${root}/parent_run
zoom_ic=${root}/zoom_l11
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-64}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
cd "${parent_run}"

if [ ! -s zoom_candidates.txt ]; then
    printf 'Missing zoom_candidates.txt\n' >&2
    exit 2
fi

python3 - zoom_candidates.txt selected_halo.txt selected_halo.args <<'PY'
import math
import pathlib
import sys

source = pathlib.Path(sys.argv[1])
report = pathlib.Path(sys.argv[2])
arguments = pathlib.Path(sys.argv[3])

rows = []
for line in source.read_text().splitlines():
    fields = line.split()
    if not fields or fields[0] == "#":
        continue
    rank = int(fields[0])
    group = int(fields[1])
    npart = int(fields[2])
    mass = float(fields[3])
    x, y, z = map(float, fields[4:7])
    rows.append(
        {"rank": rank, "group": group, "npart": npart, "mass": mass,
         "x": x, "y": y, "z": z}
    )

if not rows:
    raise SystemExit("No halo candidates were parsed")

def interior(row, margin):
    return all(margin < row[key] < 1.0 - margin for key in ("x", "y", "z"))

preferred = [
    row for row in rows
    if 2.0e13 <= row["mass"] <= 8.0e13 and interior(row, 0.08)
]
fallback = [
    row for row in rows
    if row["npart"] >= 1000 and interior(row, 0.05)
]
selected = preferred[0] if preferred else (fallback[0] if fallback else rows[0])
selection_class = "preferred_group" if preferred else (
    "interior_fallback" if fallback else "most_massive_fallback"
)

# M is in h^-1 Msun and rho_crit is in h^2 Msun Mpc^-3.  The resulting
# radius is in h^-1 Mpc.
rho_crit = 2.77536627e11
r200 = (3.0 * selected["mass"] /
        (4.0 * math.pi * 200.0 * rho_crit)) ** (1.0 / 3.0)
mask_radius_box = min(0.06, max(0.01, 4.0 * r200 / 128.0))

# Production resolution and an NFW(c=7) estimate of the Power radius.
hubble = 0.6766
zoom_particle_mass = 1.0792e10 / 8.0**3
n200_zoom = selected["mass"] / zoom_particle_mass
cell_width_l18_hinv_kpc = 128.0e3 / 2.0**18

def nfw_f(value):
    return math.log1p(value) - value / (1.0 + value)

def power_radius(kappa_target, concentration=7.0):
    lower, upper = 1.0e-7, 1.0
    for _ in range(100):
        x = 0.5 * (lower + upper)
        enclosed_fraction = nfw_f(concentration * x) / nfw_f(concentration)
        enclosed_number = max(n200_zoom * enclosed_fraction, 1.00001)
        mean_density_ratio = 200.0 * enclosed_fraction / x**3
        kappa = (
            math.sqrt(200.0) / 8.0
            * enclosed_number / math.log(enclosed_number)
            / math.sqrt(mean_density_ratio)
        )
        if kappa < kappa_target:
            lower = x
        else:
            upper = x
    return upper * r200 * 1.0e3 / hubble

report.write_text(
    "# Automatically selected parent halo\n"
    f"selection_class = {selection_class}\n"
    f"candidate_rank = {selected['rank']}\n"
    f"group = {selected['group']}\n"
    f"npart_parent = {selected['npart']}\n"
    f"approximate_mass_h-1_Msun = {selected['mass']:.8e}\n"
    f"center = {selected['x']:.10f} {selected['y']:.10f} "
    f"{selected['z']:.10f}\n"
    f"estimated_R200_h-1_Mpc = {r200:.8f}\n"
    f"mask_radius_box_units = {mask_radius_box:.10f}\n"
    "mask_extent = 4R200\n"
    f"zoom_particle_mass_h-1_Msun = {zoom_particle_mass:.8e}\n"
    f"estimated_N200_zoom = {n200_zoom:.1f}\n"
    f"level18_cell_width_h-1_kpc = {cell_width_l18_hinv_kpc:.8f}\n"
    f"level18_cell_width_physical_kpc_z0 = "
    f"{cell_width_l18_hinv_kpc / hubble:.8f}\n"
    f"estimated_Power_radius_kappa0.6_physical_kpc_c7 = "
    f"{power_radius(0.6):.8f}\n"
    f"estimated_Power_radius_kappa1.0_physical_kpc_c7 = "
    f"{power_radius(1.0):.8f}\n"
    "science_scope = resolved_group_host_core_not_kpc_subhalo_cores\n"
)
arguments.write_text(
    f"{selected['x']:.10f} {selected['y']:.10f} "
    f"{selected['z']:.10f} {mask_radius_box:.10f}\n"
)
PY

read -r xc yc zc radius < selected_halo.args
printf 'Selected halo: xc=%s yc=%s zc=%s radius=%s\n' \
       "${xc}" "${yc}" "${zc}" "${radius}"

./make_refmask.sh "${xc}" "${yc}" "${zc}" "${radius}"

point_count=$(awk 'NF >= 3 && $1 !~ /^#/' chosen_halo.part | wc -l)
if [ "${point_count}" -lt 100 ]; then
    printf 'Lagrangian mask has too few points: %s\n' "${point_count}" >&2
    exit 3
fi
printf 'Lagrangian mask points: %s\n' "${point_count}"

cd "${zoom_ic}"
if [ -d ic_zoom_l11 ]; then
    printf 'Zoom IC directory already exists. Refusing to overwrite it.\n' >&2
    exit 4
fi

./generate_zoom_ic.sh

for level in 008 009 010 011; do
    if [ ! -d "ic_zoom_l11/level_${level}" ]; then
        printf 'Missing MUSIC level_%s output\n' "${level}" >&2
        exit 5
    fi
done

cp -p music_zoom_l11.conf music_zoom_l11.used.conf
sha256sum music_zoom_l11.used.conf ../parent_run/chosen_halo.part \
          > zoom_ic_provenance.sha256
printf 'Zoom IC preparation completed at %s\n' "$(date)"
