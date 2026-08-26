#!/bin/bash
# Generate initial conditions using MUSIC
# Run from Run0 directory: bash gen_ic.sh

MUSIC=/gpfs/kjhan/Hydro/MUSIC/MUSIC2/build/MUSIC

echo "Generating IC for 128 Mpc/h, 1024^3, z=100..."
echo "Start: $(date)"

$MUSIC music.conf

echo "End: $(date)"
echo "IC generated in ic_lv10/"
