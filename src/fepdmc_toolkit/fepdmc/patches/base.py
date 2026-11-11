"""Prepare a pristine FEP-DMC pert-src (pin 05d08449) for a gfortran build with
the opt-in Rao-Blackwell window estimator.

Three groups of edits, each asserted to apply exactly once:
  1. gfortran portability : Intel
     2**31 constants, '.' component syntax, logical .eq., module wrapping of
     diagMC_gt/diagMC_JJ, HDF5 link order, MKL VSL shim. No Metropolis change.
  2. rb_window.f90 added to the build and two measurement-only hooks in
     measure_EZ_wfn (inert unless FUTUREB_RB=1).
  3. FUTUREB_SEED (opt-in) replaces the wall-clock RNG seed value, so parallel
     chains are independent and reproducible. The generator is unchanged.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

from ._anchor import locate, once, read, write

DATA = Path(__file__).resolve().parents[1] / "data"


def patch(root: str | Path, shim: str | Path | None = None) -> None:
    shim = shim or str(DATA / "mkl_vsl.f90")
    src = Path(root) / "pert-src"
    shutil.copy(shim, src / "mkl_vsl.f90")
    shutil.copy(DATA / "rb_window.f90", src / "rb_window.f90")

    p = src / "random_tool.f90"
    t = read(p)
    t = t.replace("integer, parameter ::  M = 2**31", "integer, parameter ::  M = 2147483647")
    t = t.replace(
        "integer, parameter :: A = 1103515245 , B = 12345, M = 2**31",
        "integer, parameter :: A = 1103515245 , B = 12345, M = 2147483647",
    )
    # Opt-in reproducible seed: FUTUREB_SEED replaces the wall-clock seed value.
    t = once(
        t,
        "      seed_num = 71*(mytime(7)*13+mytime(8)) + id*3499\n",
        "      seed_num = 71*(mytime(7)*13+mytime(8)) + id*3499\n"
        "      block\n"
        "         character(len=32) :: fb_env\n"
        "         integer :: fb_st, fb_seed\n"
        "         call get_environment_variable('FUTUREB_SEED', fb_env, status=fb_st)\n"
        "         if (fb_st == 0 .and. len_trim(fb_env) > 0) then\n"
        "            read(fb_env, *) fb_seed\n"
        "            seed_num = fb_seed + id*3499\n"
        "         endif\n"
        "      end block\n",
        "FUTUREB_SEED",
    )
    write(p, t)

    p = src / "diagMC.f90"
    write(
        p,
        once(
            read(p),
            "allocate(stat%update(7),stat.accept(7))",
            "allocate(stat%update(7),stat%accept(7))",
            "Intel-dot allocate",
        ),
    )
    p = src / "calc_bands.f90"
    write(p, once(read(p), "if(read_H .eq. .true.) then", "if(read_H) then", "logical .eq."))

    for name in ("diagMC_gt.f90", "diagMC_JJ.f90", "diagMC_JJt.f90", "diagMC_set.f90"):
        p = src / name
        t = read(p)
        write(p, re.sub(r"\bstat_total\.(update|accept|order|Tstat)\b", r"stat_total%\1", t))

    # RB hooks inside measure_EZ_wfn only.
    p = src / "diagMC_gt.f90"
    t = read(p)
    start = locate(t, "subroutine measure_EZ_wfn(diagram, workspace, stat)")
    end = locate(t, "end subroutine", start)
    body = t[start:end]
    body = once(
        body,
        "subroutine measure_EZ_wfn(diagram, workspace, stat)\n",
        "subroutine measure_EZ_wfn(diagram, workspace, stat)\n   use rb_window_mod, only : rb_pre, rb_post\n",
        "EZ use",
    )
    body = once(body, "   expH => workspace%expH\n", "   call rb_pre(stat)\n   expH => workspace%expH\n", "EZ pre hook")
    body = once(body, "      return \n", "      call rb_post(diagram, stat)\n      return \n", "EZ order-1 hook")
    body = once(body, "   contains \n", "   call rb_post(diagram, stat)\n   contains \n", "EZ end hook")
    t = t[:start] + body + t[end:]
    write(p, t)

    wraps = {
        "diagMC_gt.f90": ("diagmc_gt_mod", "diagmc_gt_matrix, diagmc_EZ_matrix"),
        "diagMC_JJ.f90": ("diagmc_jj_mod", "diagmc_JJ"),
    }
    uses = []
    for fname, (mod, only) in wraps.items():
        p = src / fname
        write(p, f"module {mod}\nimplicit none\ncontains\n" + read(p) + f"\nend module {mod}\n")
        uses.append(f"   use {mod}, only: {only}")
    p = src / "perturbo.f90"
    write(
        p,
        once(
            read(p),
            "   use hdf5_utils\n   implicit none",
            "   use hdf5_utils\n" + "\n".join(uses) + "\n   implicit none",
            "perturbo uses",
        ),
    )

    p = src / "makefile"
    t = read(p)
    t = once(t, "diagMC_updates.f90 \n", "diagMC_updates.f90 \\\nrb_window.f90 \n", "makefile rb object")
    t = once(t, "$(QELIBS) -lstdc++", "$(QELIBS) $(HDF5_LIB) -lstdc++", "HDF5 link order")
    write(p, t)

    shutil.copy(DATA / "make.sys", Path(root) / "make.sys")
