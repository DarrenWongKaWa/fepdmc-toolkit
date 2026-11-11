"""Add the grouped-measure chain (bchain.f90) to a tree already prepared by
patches/base.py. Hooks go only into diagmc_EZ_matrix, each asserted once:
around update_drive (second-stage correction, mode refresh) and around
measure_EZ_wfn (grouped estimator). Inert unless FUTUREB_BCHAIN=1.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from ._anchor import locate, once, read, write

DATA = Path(__file__).resolve().parents[1] / "data"


def patch(root: str | Path) -> None:
    src = Path(root) / "pert-src"
    shutil.copy(DATA / "bchain.f90", src / "bchain.f90")
    p = src / "diagMC_gt.f90"
    t = read(p)
    start = locate(t, "subroutine diagmc_EZ_matrix()")
    end = locate(t, "end subroutine", start)
    body = t[start:end]
    body = once(
        body,
        "   use multiphonon_update_matrix, only : update_drive\n",
        "   use multiphonon_update_matrix, only : update_drive\n"
        "   use bchain_mod, only : bchain_pre_update, bchain_post_update, "
        "bchain_pre_measure, bchain_post_measure\n",
        "EZ use",
    )
    body = once(
        body,
        "         call update_drive(diagrams(i_omp),stat(i_omp),update)\n",
        "         call bchain_pre_update(diagrams(i_omp), stat(i_omp))\n"
        "         call update_drive(diagrams(i_omp),stat(i_omp),update)\n"
        "         call bchain_post_update(diagrams(i_omp), stat(i_omp))\n",
        "EZ update hooks",
    )
    body = once(
        body,
        "            call measure_EZ_wfn(diagrams(i_omp),workspaces(i_omp),stat(i_omp))\n",
        "            call bchain_pre_measure(stat(i_omp))\n"
        "            call measure_EZ_wfn(diagrams(i_omp),workspaces(i_omp),stat(i_omp))\n"
        "            call bchain_post_measure(diagrams(i_omp), stat(i_omp))\n",
        "EZ measure hooks",
    )
    write(p, t[:start] + body + t[end:])
    p = src / "makefile"
    write(p, once(read(p), "rb_window.f90 \n", "rb_window.f90 \\\nbchain.f90 \n", "makefile bchain"))
