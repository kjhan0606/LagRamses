# Driver disposition: static Fe spectral connection

Fable's single read-only review returned REVISE. It identified the real
CHIMES/layout conflict, and correctly warned that a broad 1--5.6 eV maxent
packet violates the retained 4 eV absorption limit even for a low mean.
Q-GOAL was conditional on a useful runnable source; Q-LEAN judged the
implementation compact. The raw review is retained separately. No second
review/end-audit is required by the operator's selective-audit policy.

Adopt finding1: CHIMES-free Fe is now ONLY static seed mass: cooling=none,
zero C/silicate/Fe condensation, no growth/sputtering/size transfer/SN shocks,
no Fe kinetics or relative motion. Two Fe fields follow the reserved dust
window; existing CHIMES offsets remain unchanged. Make's DUST_IRON profile
is NVAR32 without CHIMES and the old NVAR189 with CHIMES (NENER1/virial).
Both mkrun and generic namelist/GUI validation follow these prerequisites.

Choose the first explicit option for finding2: retain the approved <=4 eV
and <=300 K bounds. The suggestion to assign all hard-photon absorption to
grain heat is NOT adopted merely to make the run pass. The new mode is
therefore a bounded static/sub-eV comparison, NOT general stellar/AGN Fe
spectral qualification and NOT completion of medium item8. An endpoint
packet in group2 can remain below4 eV, but a broad mean2.366 eV packet must
reject its actual >4 eV tail; no clipping, renormalization or tail tolerance.
The live test uses an explicitly synthetic group1 source, not a truncated
SED misrepresented as physical. Native tests cover nonzero Fe opacity,
actual energy, Fe-only/mixed cells, opaque limit and hard-tail rollback.

Optional items adopted: bitwise zero-Fe equality against the original
four-bin C ABI, and distinct Fe-limit status10 instead of secondary status9.
No shared-builder refactor/new harness. The actual128-versus2048 Mie check
is retained as one numerical measurement; missing spin response, frozen
dielectric and shared material temperature remain explicit limitations.

This bounded extension exercises/reuses existing native physics without
removing scientific admission limits. A genuinely general Fe stellar/AGN
mode still needs a physically specified photoelectron/charge/energy closure;
the static connection does not resolve that scientific source gap.
