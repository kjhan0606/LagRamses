#!/bin/bash
# [RESIZABLE] Expand a Linux CPU-list mask into a canonical CSV list.
# Usage: expand_cpuset.sh <cpu-list> <expected-count>

set -euo pipefail

if [ "$#" -ne 2 ] || [[ ! $2 =~ ^[1-9][0-9]*$ ]]; then
  echo "usage: $0 <cpu-list> <expected-count>" >&2
  exit 64
fi

mask=$1
expected=$2

awk -v mask="$mask" -v expected="$expected" '
  BEGIN {
    nchunk=split(mask,chunk,",")
    count=0
    previous=-1
    if(nchunk<1 || mask=="") exit 65

    for(i=1;i<=nchunk;i++) {
      if(chunk[i] ~ /^[0-9]+$/) {
        first=chunk[i]+0
        last=first
      } else {
        nrange=split(chunk[i],range,"-")
        if(nrange!=2 || range[1] !~ /^[0-9]+$/ || \
           range[2] !~ /^[0-9]+$/) exit 66
        first=range[1]+0
        last=range[2]+0
        if(first>last) exit 67
      }

      for(cpu=first;cpu<=last;cpu++) {
        # Canonical input must be strictly increasing.  This also rejects
        # duplicates and overlapping ranges instead of silently remapping.
        if(seen[cpu]++ || cpu<=previous) exit 68
        ids[++count]=cpu
        previous=cpu
      }
    }

    if(count!=expected) exit 69
    for(i=1;i<=count;i++)
      printf "%s%d", (i==1 ? "" : ","), ids[i]
    printf "\n"
  }
' </dev/null
