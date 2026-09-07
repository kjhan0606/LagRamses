"""Offline He-accretion reference; NOT an N100 source or binary evolution solver.

Wang, Podsiadlowski & Han 2017, equations 1--2 (arXiv:1708.07067):
solar-composition, hot, nonrotating CO WDs. KH04 (astro-ph/0407632),
equations 1--6: conditional flash retention, no unmodelled Roche-lobe loss.
Mass interpolation of KH04 node fits is a declared comparison approximation.
No extrapolation, artificial wind cap, or conversion of off-centre ignition
to a prescribed NS/OSi fate is made here. A local regime is not an SN event.
"""
from __future__ import annotations

import bisect
import math

MODEL_ID = 'W17_stable_KH04_flash_conditional_v1'
MASSES = (.7, .8, .9, 1., 1.1, 1.2, 1.3, 1.35)
LOWER_LOG_RATE = (-7.4, -6.5, -6.88, -6.92, -7.06, -7.06, -7.35, -7.4)


def _kh_node(index, log_rate):
    if log_rate <= LOWER_LOG_RATE[index]:
        return None
    mass = MASSES[index]
    if mass == .7:
        raw = 1.  # KH04 text: only if the expanded envelope fits its Roche lobe.
    elif mass == .8:
        raw = 1. if log_rate >= -6.34 else -.35*(log_rate+6.1)**2+1.02
    elif mass == .9:
        raw = 1. if log_rate >= -6.05 else -.35*(log_rate+5.6)**2+1.07
    elif mass == 1.:
        raw = 1. if log_rate >= -5.93 else -.35*(log_rate+5.6)**2+1.01
    elif mass in (1.1, 1.2):
        if log_rate < -5.95:
            raw = .54*log_rate+4.16
        elif log_rate < -5.76:
            raw = -.54*(log_rate+5.6)**2+1.01
        else:
            raw = 1.
    elif mass == 1.3:
        raw = 1. if log_rate >= -5.83 else -.175*(log_rate+5.35)**2+1.03
    else:
        raw = 1. if log_rate >= -6.05 else -.115*(log_rate+5.7)**2+1.01
    # Printed polynomial joins can slightly exceed one. This explicit bound
    # only limits KH04's positive retention fits, not general erosion models.
    return min(raw, 1.)


def kh04_conditional_efficiency(mass, mdot):
    """Return None outside common node support; never substitute zero retention."""
    if not math.isfinite(mass) or not math.isfinite(mdot) or mdot <= 0:
        raise ValueError('finite mass and positive finite rate required')
    if not MASSES[0] <= mass <= MASSES[-1]:
        return None
    x = math.log10(mdot)
    # The branch formulae also describe the high-retention side; the caller
    # must decide steady burning/expansion using W17, not KH04 alone.
    upper = bisect.bisect_left(MASSES, mass)
    if MASSES[upper] == mass:
        return _kh_node(upper, x)
    lower = upper-1
    a,b = _kh_node(lower,x),_kh_node(upper,x)
    if a is None or b is None:
        return None
    f = (mass-MASSES[lower])/(MASSES[upper]-MASSES[lower])
    return (1-f)*a+f*b


def helium_regime(mass, mdot):
    """Classify the incident He rate in Msun/yr under the reference assumptions.

    Does NOT validate a supplied donor-loss/net-WD-growth rate as incident He.
    Does NOT imply N100 compatibility of the evolving core or event mass.
    """
    if not math.isfinite(mass) or not math.isfinite(mdot) or mdot <= 0:
        raise ValueError('finite mass and positive finite rate required')
    result = dict(model_id=MODEL_ID, regime='outside_W17_grid', eta_conditional=None,
                  stable_min_msun_yr=None, expansion_min_msun_yr=None,
                  offcentre_risk_reference=bool(mdot > 2.05e-6),
                  N100_event_qualified=False)
    if not .6 <= mass <= 1.35 or not 1e-8 <= mdot <= 1e-5:
        return result
    stable = 1.46e-6*(-mass**3+3.45*mass**2-2.60*mass+.85)
    expansion = 2.17e-6*(mass**2+.82*mass-.38)
    result.update(stable_min_msun_yr=stable, expansion_min_msun_yr=expansion)
    if mdot > expansion:
        result['regime'] = 'expansion_requires_coupled_wind_or_RLO_model'
    elif mdot >= stable:
        result['eta_conditional'] = 1.
        result['regime'] = ('stable_offcentre_risk_not_N100' if mdot > 2.05e-6
                            else 'stable_central_growth_candidate_only')
    else:
        eta = kh04_conditional_efficiency(mass,mdot)
        result['eta_conditional'] = eta
        result['regime'] = ('flash_outside_KH04_support' if eta is None
                            else 'flash_retention_requires_Roche_geometry')
    return result
