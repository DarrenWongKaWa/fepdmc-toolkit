"""fepdmc-toolkit: tools for first-principles electron-phonon diagrammatic Monte Carlo.

Subpackages:

  fepdmc            source patches for an upstream FEP-DMC tree (opt-in fixes, exact
                    moves, diagnostics), a Docker run driver, the pooled ratio estimator,
                    the end-to-end exactness check and analyses of run directories;
                    command line ``fepdmc-toolkit``
  p1_da             algebra of exact-corrected delayed acceptance (P1)
  formation_energy  formation energy Q = ratio - E_bare
  stats             jackknife, HAC and bootstrap helpers for correlated Monte Carlo
  event_trace       parsers for the counters printed by the patched native code
  fortran           Fortran sources of the benchmarked P1 build (reference, GPL-3)
  diagram_compiler  compiler from fixed-order diagram families to shared DAGs, a
                    delayed-acceptance kernel and generated Fortran
  scheme_d          reference IR, interpreter and exact oracles used to validate the compiler
  r1_fixed_order    recorded result of the fixed-order R1 evaluator benchmark
  r1_grouped        finite signed toy model with two grouped measures and a dense oracle

The package does not include FEP-DMC, Perturbo or Quantum ESPRESSO.
"""

import importlib

__version__ = "1.0.0"

_LAZY = {
    "E_BARE_LIF_METHOD0": "formation_energy",
    "formation_energy_eV": "formation_energy",
    "delayed_acceptance_prob": "p1_da",
    "native_re_ratio": "p1_da",
}

__all__ = list(_LAZY)


def __getattr__(name):
    if name in _LAZY:
        return getattr(importlib.import_module(f".{_LAZY[name]}", __name__), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
