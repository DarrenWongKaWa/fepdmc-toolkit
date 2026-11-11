"""Opt-in fix for the reference trace in FEP-DMC's remove_external_ph.

When the diagram has internal vertices, remove_external_ph builds the trace of
the configuration *without* the external pair (mat_old) with the wrap segment
propagating at v_head%ekout. At that point the pair is still present, so the head
segment carries momentum k - q. The pair-less wrap segment carries k, whose
energies are vn1%ekout (= vL%ekin). The removal acceptance therefore compares the
with-pair trace with a reference on the wrong band energies. add_external_ph's
identical line is correct (there the pair does not exist yet), and so is the
empty-diagram branch of the removal, which already uses vn1%ekout.

With FUTUREB_EXTRMFIX=1 the reference wrap segment uses vn1%ekout; otherwise the
code is unchanged.
"""

from __future__ import annotations

from pathlib import Path

from ._anchor import locate, once, read, write


def patch(root: str | Path) -> None:
    p = Path(root) / "pert-src" / "diagMC_JJ_updates.f90"
    t = read(p)
    mod = locate(t, "module multiphonon_update_matrix")
    start = locate(t, "subroutine remove_external_ph(diagram, stat)", mod)
    end = locate(t, "end subroutine", start)
    body = t[start:end]
    body = once(
        body,
        "      type(vertex),pointer :: vL,vR, vn1, vn2 \n",
        "      type(vertex),pointer :: vL,vR, vn1, vn2 \n"
        "      logical, save :: rmfix_checked = .false., rmfix = .false.\n"
        "      character(len=8) :: rmenv\n"
        "      integer :: rmst\n",
        "declarations",
    )
    body = once(
        body,
        "         !old : vR -> v_tail/v_head -> vL\n         Elist = v_head%ekout - minval(v_head%ekout)\n",
        "         !old : vR -> v_tail/v_head -> vL\n"
        "         if (.not. rmfix_checked) then\n"
        "            call get_environment_variable('FUTUREB_EXTRMFIX', rmenv, status=rmst)\n"
        "            rmfix = (rmst == 0 .and. trim(adjustl(rmenv)) == '1')\n"
        "            rmfix_checked = .true.\n"
        "         end if\n"
        "         if (rmfix) then\n"
        "            Elist = vn1%ekout - minval(vn1%ekout)\n"
        "         else\n"
        "            Elist = v_head%ekout - minval(v_head%ekout)\n"
        "         end if\n",
        "reference wrap segment",
    )
    write(p, t[:start] + body + t[end:])
