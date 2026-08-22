#!/usr/bin/env python3
"""Focused fail-closed tests for CUDA/nGR numerical comparison helpers."""

from __future__ import annotations

import math
import unittest

import compare_cuda_ndgp_outputs as comparator


class ParticleFiniteTest(unittest.TestCase):
    @staticmethod
    def particles() -> dict[int, comparator.Particle]:
        result: dict[int, comparator.Particle] = {}
        for identity in comparator.pinned_particle_ids():
            result[identity] = comparator.Particle(
                (0.01, 0.02, 0.03),
                (0.1, 0.2, 0.3),
                1.0,
                5 if identity <= 32**3 else 6,
                1,
                -0.5,
            )
        return result

    def test_particle_nan_is_rejected(self) -> None:
        cpu = self.particles()
        gpu = dict(cpu)
        identity = min(gpu)
        original = gpu[identity]
        gpu[identity] = comparator.Particle(
            (math.nan, original.position[1], original.position[2]),
            original.velocity,
            original.mass,
            original.level,
            original.particle_type,
            original.potential,
        )
        with self.assertRaises(comparator.FormatError):
            comparator.compare_particles(cpu, gpu)


if __name__ == "__main__":
    unittest.main()
