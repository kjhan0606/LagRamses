#!/bin/bash
# [RESIZABLE] Validate MPI rank/localid/single-CPU rows against an expected map.
# Usage: validate_affinity_map.sh <expected-cpu-csv> <map-file> <rank-count>

set -euo pipefail

if [ "$#" -ne 3 ] || [[ ! $3 =~ ^[1-9][0-9]*$ ]] || [ ! -r "$2" ]; then
  echo "usage: $0 <expected-cpu-csv> <map-file> <rank-count>" >&2
  exit 64
fi

expected=$1
map_file=$2
nproc=$3

awk -v expected="$expected" -v nproc="$nproc" '
  BEGIN {
    nexpected=split(expected,cpu,",")
    if(nexpected!=nproc) bad=1
    for(i=1;i<=nexpected;i++) {
      if(cpu[i] !~ /^[0-9]+$/ || expected_cpu[cpu[i]+0]++) bad=1
    }
  }
  NF!=3 || $1 !~ /^[0-9]+$/ || $2 !~ /^[0-9]+$/ || \
  $3 !~ /^[0-9]+$/ {
    bad=1
    next
  }
  {
    if($1<0 || $1>=nproc || $2<0 || $2>=nproc) {
      bad=1
      next
    }
    rank[$1]++
    localid[$2]++
    mask[$3+0]++
    if(($3+0)!=(cpu[$2+1]+0)) bad=1
  }
  END {
    if(NR!=nproc) bad=1
    for(i=0;i<nproc;i++) {
      if(rank[i]!=1) bad=1
      if(localid[i]!=1) bad=1
    }
    nmask=0
    for(m in mask) nmask++
    if(nmask!=nproc) bad=1
    exit(bad ? 1 : 0)
  }
' "$map_file"
