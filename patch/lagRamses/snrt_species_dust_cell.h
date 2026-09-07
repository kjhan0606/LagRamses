#pragma once
// One scientific cell operator shared by CUDA and OpenMP. Group order,
// inventory guard band and photon/dust ledgers are intentionally identical.
#include <cfloat>
#include <cmath>
#ifdef __CUDACC__
#define SNRT_CELL_HD __host__ __device__
#else
#define SNRT_CELL_HD
#endif
SNRT_CELL_HD inline void snrt_cap_species_dust_cell(
    float *state, float *absorbed_direction,
    const float *optical_depth_species, const float *optical_depth_dust,
    float *available_species, float *absorbed_hhe_species,
    float *absorbed_dust_group, float *returned_group, float *raw_group,
    float *absorbed_group, float *absorbed_total,
    int nowned, int nwork, int ndirection, int ngroup, int cell) {
  if (cell >= nowned) return;

  const long long group_count = static_cast<long long>(nowned) * ngroup;
  float available[3];
  for (int species = 0; species < 3; ++species) {
    available[species] = available_species[species * nowned + cell];
  }
  for (int group = 0; group < ngroup; ++group) {
    const long long group_base = static_cast<long long>(group) * nwork * ndirection;
    const long long output_index = static_cast<long long>(group) * nowned + cell;
    const long long hhe_base = static_cast<long long>(group) * nowned + cell;
    float opacity[3];
    float hhe_tau = 0.0f;
    for (int species = 0; species < 3; ++species) {
      opacity[species] = optical_depth_species[species * group_count + hhe_base];
      hhe_tau += opacity[species];
      absorbed_hhe_species[species * group_count + hhe_base] = 0.0f;
    }
    const float dust_tau = optical_depth_dust[output_index];
    const float component_tau = hhe_tau + dust_tau;
    float raw_absorbed = 0.0f;
    for (int idir = 0; idir < ndirection; ++idir) {
      raw_absorbed += fmaxf(0.0f,
          absorbed_direction[group_base + static_cast<long long>(idir) * nwork + cell]);
    }
    raw_group[output_index] = raw_absorbed;
    if (raw_absorbed <= 0.0f) {
      absorbed_dust_group[output_index] = 0.0f;
      returned_group[output_index] = 0.0f;
      absorbed_group[output_index] = 0.0f;
      continue;
    }

    float hhe_target_full = 0.0f;
    if (dust_tau == 0.0f) {
      // Exact legacy arithmetic in the zero-dust limit.
      hhe_target_full = raw_absorbed;
    } else if (hhe_tau > 0.0f && component_tau > 0.0f) {
      hhe_target_full = fminf(raw_absorbed,
          raw_absorbed * hhe_tau / component_tau);
    }
    hhe_target_full = fmaxf(0.0f, fminf(raw_absorbed, hhe_target_full));

    float eligible_inventory = 0.0f;
    for (int species = 0; species < 3; ++species) {
      if (opacity[species] > 0.0f) eligible_inventory += available[species];
    }
    float inventory_scale = fmaxf(raw_absorbed, eligible_inventory);
    for (int species = 0; species < 3; ++species) {
      inventory_scale = fmaxf(inventory_scale, fabsf(available[species]));
    }
    const float remainder_tolerance =
        256.0f * FLT_EPSILON * fmaxf(inventory_scale, FLT_MIN);
    // Match the legacy guard band.  Its slice is explicitly returned below,
    // not treated as physical H/He excess eligible for dust transfer.
    const float target_absorbed = fminf(hhe_target_full, eligible_inventory * 0.99995f);
    float assigned[3] = {0.0f, 0.0f, 0.0f};
    float remaining = target_absorbed;
    bool active[3] = {opacity[0] > 0.0f && available[0] > 0.0f,
                      opacity[1] > 0.0f && available[1] > 0.0f,
                      opacity[2] > 0.0f && available[2] > 0.0f};
    float active_weight = 0.0f;
    for (int species = 0; species < 3; ++species) {
      if (active[species]) active_weight += opacity[species];
    }
    for (int pass = 0; pass < 3 && remaining > remainder_tolerance &&
         active_weight > 0.0f; ++pass) {
      const float pass_remaining = remaining;
      const float pass_weight = active_weight;
      bool saturated_any = false;
      bool saturated[3] = {false, false, false};
      for (int species = 0; species < 3; ++species) {
        if (!active[species]) continue;
        const float headroom = fmaxf(0.0f, available[species] - assigned[species]);
        if (headroom <= remainder_tolerance) {
          saturated[species] = true;
          saturated_any = true;
          continue;
        }
        const float requested = pass_remaining * opacity[species] / pass_weight;
        const float addition = fminf(requested, headroom);
        assigned[species] += addition;
        if (addition + remainder_tolerance >= headroom) {
          saturated[species] = true;
          saturated_any = true;
        }
      }
      remaining = target_absorbed - (assigned[0] + assigned[1] + assigned[2]);
      if (remaining < 0.0f) remaining = 0.0f;
      if (!saturated_any) {
        remaining = 0.0f;
        break;
      }
      for (int species = 0; species < 3; ++species) {
        if (saturated[species]) {
          active[species] = false;
          active_weight -= opacity[species];
        }
      }
    }
    if (remaining > remainder_tolerance) {
      for (int species = 0; species < 3 && remaining > 0.0f; ++species) {
        if (!active[species]) continue;
        const float headroom = fmaxf(0.0f, available[species] - assigned[species]);
        const float addition = fminf(headroom, remaining);
        assigned[species] += addition;
        remaining -= addition;
      }
    }

    const float assigned_hhe = assigned[0] + assigned[1] + assigned[2];
    const float physical_excess = fmaxf(0.0f, hhe_target_full - eligible_inventory);
    const float guard_return = fmaxf(0.0f,
        hhe_target_full - target_absorbed - physical_excess);
    const float dust_direct = fmaxf(0.0f, raw_absorbed - hhe_target_full);
    const float dust_fraction = -expm1f(-dust_tau);
    const float dust_transfer = physical_excess * dust_fraction;
    const float hhe_residual = fmaxf(0.0f, target_absorbed - assigned_hhe);
    // Keep the guard-band and the unabsorbed part of the finite H/He excess
    // explicit in the returned ledger.  Defining the assigned total as
    // raw-returned makes the group closure exact in FP32; the dust ledger is
    // then the non-H/He remainder of that assigned total.
    float returned = 0.0f;
    float assigned_total = 0.0f;
    float assigned_dust = 0.0f;
    if (dust_tau == 0.0f) {
      // Keep the old cap and directional arithmetic exactly in the zero-dust
      // limit.  The returned ledger is the old unassigned remainder.
      assigned_total = assigned_hhe;
      returned = fmaxf(0.0f, raw_absorbed - assigned_total);
    } else {
      const float returned_model = guard_return + hhe_residual +
          fmaxf(0.0f, physical_excess - dust_transfer);
      returned = fminf(raw_absorbed, fmaxf(0.0f, returned_model));
      assigned_total = fmaxf(0.0f, raw_absorbed - returned);
      assigned_dust = fmaxf(0.0f, assigned_total - assigned_hhe);
      if (assigned_hhe > assigned_total) {
        assigned_total = assigned_hhe;
        assigned_dust = 0.0f;
      }
    }
    for (int species = 0; species < 3; ++species) {
      absorbed_hhe_species[species * group_count + hhe_base] = assigned[species];
      available[species] = fmaxf(0.0f, available[species] - assigned[species]);
    }
    const float cap = assigned_total > 0.0f ? assigned_total / raw_absorbed : 0.0f;
    for (int idir = 0; idir < ndirection; ++idir) {
      const long long index = group_base + static_cast<long long>(idir) * nwork + cell;
      const float removed = fmaxf(0.0f, absorbed_direction[index]);
      const float limited_removed = removed * cap;
      state[index] += removed - limited_removed;
      absorbed_direction[index] = limited_removed;
    }
    // Match the legacy per-group reduction order.  This keeps the zero-dust
    // ABI comparison meaningful while the raw/returned ledgers remain
    // independently available to the DUST-8 host partition.
    float assigned_group = 0.0f;
    for (int idir = 0; idir < ndirection; ++idir) {
      assigned_group += absorbed_direction[group_base +
          static_cast<long long>(idir) * nwork + cell];
    }
    absorbed_dust_group[output_index] = assigned_dust;
    returned_group[output_index] = fmaxf(0.0f, raw_absorbed - assigned_group);
    absorbed_group[output_index] = assigned_group;
  }

  // Reproduce snrt_reduce_multigroup_kernel's 128-lane reduction order for
  // the cell-total output.  This is intentionally separate from the ledger
  // arithmetic so the zero-dust regression can compare the legacy scalar
  // absorption bit-for-bit.
  float partial[128];
  for (int lane = 0; lane < 128; ++lane) {
    float sum = 0.0f;
    for (int group = 0; group < ngroup; ++group) {
      const long long group_base = static_cast<long long>(group) * nwork * ndirection;
      for (int idir = lane; idir < ndirection; idir += 128) {
        sum += absorbed_direction[group_base + static_cast<long long>(idir) * nwork + cell];
      }
    }
    partial[lane] = sum;
  }
  for (int stride = 64; stride > 0; stride /= 2) {
    for (int lane = 0; lane < stride; ++lane) {
      partial[lane] += partial[lane + stride];
    }
  }
  absorbed_total[cell] = partial[0];

  for (int species = 0; species < 3; ++species) {
    available_species[species * nowned + cell] = available[species];
  }
}
#undef SNRT_CELL_HD
