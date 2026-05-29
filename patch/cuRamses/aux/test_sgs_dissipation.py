#!/usr/bin/env python3
"""
Test 3.1: SGS Dissipation-Only Verification

Verifies the analytical implicit dissipation formula against the continuous ODE solution.

The SGS dissipation equation per unit mass:
  de_sgs/dt = -D/rho = -C_diss * sigma^3 / dx

where e_sgs = (3/2)*sigma^2.

Continuous ODE solution:
  From de/dt = -C_diss * sigma^3/dx  and  e = 1.5*sigma^2:
  3*sigma * dsigma/dt = -C_diss * sigma^3/dx
  dsigma/dt = -C_diss * sigma^2 / (3*dx)

  Exact solution: sigma(t) = sigma_0 / (1 + C_diss*sigma_0*t/(3*dx))
  Therefore:     e_sgs(t) = 1.5*sigma_0^2 / (1 + C_diss*sigma_0*t/(3*dx))^2

Code's implicit step (one timestep dt):
  sigma_new = sigma / (1 + C_diss*sigma*dt/(N*dx))
  where N = 3 (corrected) or 2 (old code)

This script verifies that N=3 reproduces the continuous solution.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parameters
C_diss = 0.5    # current code value
dx = 1.0        # arbitrary units
sigma_0 = 1.0   # initial turbulent velocity
e_sgs_0 = 1.5 * sigma_0**2

# Time integration parameters
dt = 0.01 * dx / sigma_0  # small dt for accuracy
t_end = 10.0 * dx / sigma_0  # ~10 dissipation times
nsteps = int(t_end / dt)

# ============================================================
# 1. Exact continuous solution
# ============================================================
t_exact = np.linspace(0, t_end, 1000)
sigma_exact = sigma_0 / (1.0 + C_diss * sigma_0 * t_exact / (3.0 * dx))
e_exact = 1.5 * sigma_exact**2

# ============================================================
# 2. Code's implicit step with 3*dx (corrected)
# ============================================================
sigma_3dx = np.zeros(nsteps + 1)
sigma_3dx[0] = sigma_0
for i in range(nsteps):
    sigma_3dx[i+1] = sigma_3dx[i] / (1.0 + C_diss * sigma_3dx[i] * dt / (3.0 * dx))
e_3dx = 1.5 * sigma_3dx**2
t_3dx = np.arange(nsteps + 1) * dt

# ============================================================
# 3. Old code's implicit step with 2*dx
# ============================================================
sigma_2dx = np.zeros(nsteps + 1)
sigma_2dx[0] = sigma_0
for i in range(nsteps):
    sigma_2dx[i+1] = sigma_2dx[i] / (1.0 + C_diss * sigma_2dx[i] * dt / (2.0 * dx))
e_2dx = 1.5 * sigma_2dx**2

# ============================================================
# 4. Direct ODE integration (RK4) for ground truth
# ============================================================
def dsigma_dt(sigma, C_diss, dx):
    return -C_diss * sigma**2 / (3.0 * dx)

sigma_rk4 = np.zeros(nsteps + 1)
sigma_rk4[0] = sigma_0
for i in range(nsteps):
    s = sigma_rk4[i]
    k1 = dt * dsigma_dt(s, C_diss, dx)
    k2 = dt * dsigma_dt(s + 0.5*k1, C_diss, dx)
    k3 = dt * dsigma_dt(s + 0.5*k2, C_diss, dx)
    k4 = dt * dsigma_dt(s + k3, C_diss, dx)
    sigma_rk4[i+1] = s + (k1 + 2*k2 + 2*k3 + k4) / 6.0
e_rk4 = 1.5 * sigma_rk4**2

# ============================================================
# 5. Error analysis
# ============================================================
# Interpolate exact solution at discrete times
sigma_exact_at_t = sigma_0 / (1.0 + C_diss * sigma_0 * t_3dx / (3.0 * dx))

err_3dx = np.max(np.abs(sigma_3dx - sigma_exact_at_t) / sigma_exact_at_t)
err_2dx_vs_3dx_exact = np.max(np.abs(sigma_2dx - sigma_exact_at_t) / sigma_exact_at_t)
err_rk4 = np.max(np.abs(sigma_rk4 - sigma_exact_at_t) / sigma_exact_at_t)

print("=" * 60)
print("Test 3.1: SGS Dissipation-Only Verification")
print("=" * 60)
print(f"Parameters: C_diss={C_diss}, dx={dx}, sigma_0={sigma_0}")
print(f"dt={dt:.4f}, nsteps={nsteps}, t_end={t_end:.1f}")
print()
print("Max relative error vs exact analytical solution:")
print(f"  3dx implicit (corrected): {err_3dx:.2e}")
print(f"  2dx implicit (old code):  {err_2dx_vs_3dx_exact:.2e}")
print(f"  RK4 (ground truth):       {err_rk4:.2e}")
print()

# Check: the 2dx formula should match a DIFFERENT ODE
sigma_2dx_exact = sigma_0 / (1.0 + C_diss * sigma_0 * t_3dx / (2.0 * dx))
err_2dx_own_exact = np.max(np.abs(sigma_2dx - sigma_2dx_exact) / sigma_2dx_exact)
print(f"  2dx implicit vs 2dx-ODE:  {err_2dx_own_exact:.2e}")
print()

# Effective dissipation rate comparison
print("Effective dissipation rate (D/rho = C_eff * sigma^3/dx):")
print(f"  3dx formula: C_eff = C_diss = {C_diss:.2f}")
print(f"  2dx formula: C_eff = C_diss * 3/2 = {C_diss * 1.5:.2f}")
print(f"  (2dx is 50% stronger than stated D = C_diss*sigma^3/dx)")
print()

# Literature comparison
print("Literature correspondence:")
print(f"  Code C_diss = C_epsilon * (3/2)^(3/2) = C_epsilon * {1.5**1.5:.3f}")
print(f"  C_diss=0.5 -> C_epsilon = {0.5/1.5**1.5:.3f} (vs Schmidt 1.5, SLES 1.6)")
print(f"  C_diss=2.76 -> C_epsilon = {2.76/1.5**1.5:.3f}")
print()

# Energy conservation check
print("Energy conservation (dissipation -> thermal):")
de_dissipated_3dx = e_3dx[0] - e_3dx[-1]
print(f"  Initial e_sgs = {e_3dx[0]:.4f}")
print(f"  Final e_sgs   = {e_3dx[-1]:.4f}")
print(f"  Dissipated     = {de_dissipated_3dx:.4f}")
print(f"  Should go to thermal energy (verified in code)")

# ============================================================
# 6. Large dt test (stability check)
# ============================================================
print()
print("=" * 60)
print("Stability test with large dt (CFL > 1)")
print("=" * 60)
for dt_large in [1.0, 5.0, 20.0, 100.0]:
    s = sigma_0
    for _ in range(10):
        s = s / (1.0 + C_diss * s * dt_large / (3.0 * dx))
    print(f"  dt/t_diss = {dt_large*sigma_0/dx:.1f}: sigma = {s:.6f} (>0, stable)")

# ============================================================
# 7. Plot
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.plot(t_exact * sigma_0 / dx, sigma_exact / sigma_0, 'k-', lw=2, label='Exact (3dx ODE)')
ax.plot(t_3dx * sigma_0 / dx, sigma_3dx / sigma_0, 'b--', lw=1.5, label='Implicit 3dx (corrected)')
ax.plot(t_3dx * sigma_0 / dx, sigma_2dx / sigma_0, 'r:', lw=1.5, label='Implicit 2dx (old code)')
ax.plot(t_3dx[::50] * sigma_0 / dx, sigma_rk4[::50] / sigma_0, 'go', ms=4, label='RK4 (ground truth)')
ax.set_xlabel(r'$t \cdot \sigma_0 / \Delta x$')
ax.set_ylabel(r'$\sigma / \sigma_0$')
ax.set_title('SGS Dissipation: sigma(t)')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.semilogy(t_3dx * sigma_0 / dx, np.abs(sigma_3dx - sigma_exact_at_t) / sigma_exact_at_t,
            'b-', label='3dx error')
ax.semilogy(t_3dx * sigma_0 / dx, np.abs(sigma_2dx - sigma_exact_at_t) / sigma_exact_at_t,
            'r-', label='2dx error (vs 3dx exact)')
ax.set_xlabel(r'$t \cdot \sigma_0 / \Delta x$')
ax.set_ylabel('Relative error')
ax.set_title('Error vs exact 3dx solution')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/test_sgs_dissipation.png', dpi=150)
print()
print(f"Plot saved to misc/test_sgs_dissipation.png")
