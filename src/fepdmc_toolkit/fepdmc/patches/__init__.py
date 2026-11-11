"""Source patches for a pristine FEP-DMC tree (upstream pin in fepdmc_toolkit.fepdmc.FEP_DMC_PIN).

Each module exposes ``patch(root)`` and asserts every anchor it edits exactly once,
so a patch either applies cleanly to the pinned tree or raises :class:`PatchError`.

  base      gfortran portability, MKL VSL shim, measurement hooks, FUTUREB_SEED
  bchain    toolkit moves and diagnostics (bchain.f90; inert unless switched on)
  wqfix     fix: stale phonon frequency in add_external_ph   (FUTUREB_WQFIX=1)
  extfix    fix: reverse-proposal support in remove_external_ph (FUTUREB_EXTFIX=1)
  extrmfix  fix: reference trace in remove_external_ph       (FUTUREB_EXTRMFIX=1)

The three fixes are opt-in at run time: without the switch the code path is the
upstream one.

``apply`` edits the tree in place, one patch at a time. Use
``fepdmc_toolkit.fepdmc.prepare.prepare_tree`` (command ``fepdmc-toolkit prepare``)
to apply a whole profile: it tries the patches on a copy first and writes the tree
only if all of them apply.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path

from ._anchor import PatchError

__all__ = ["ORDER", "PROFILES", "PatchError", "apply"]

ORDER = ("base", "bchain", "wqfix", "extfix", "extrmfix")

PROFILES = {
    "fixes": ("base", "wqfix", "extfix", "extrmfix"),
    "research": ORDER,
}


def apply(root: str | Path, name: str) -> None:
    """Apply one named patch to the perturbo-fep-dmc directory ``root``, in place."""
    if name not in ORDER:
        raise ValueError(f"unknown patch {name!r}; expected one of {ORDER}")
    import_module(f"{__name__}.{name}").patch(root)
