"""FEP-DMC patches, run tools and estimators.

Works on an upstream FEP-DMC checkout at ``FEP_DMC_PIN`` (Luo, Park, Bernardi). It
does not ship FEP-DMC itself; patches are applied to a user-provided tree and built
in the user's QE 6.5 toolchain.

  prepare    one-step patching of a pristine tree (profiles "fixes", "research");
             the tree is written only if every patch applies
  patches    the individual patches; the three fixes are opt-in at run time
             (FUTUREB_WQFIX, FUTUREB_EXTFIX, FUTUREB_EXTRMFIX). WQFIX is the defect
             that the C0 adapter (integration/c0) fixes unconditionally, made switchable.
  runner     Docker driver for diagmc-EZ chains with the toolkit switches
  pooled     pooled ratio estimator over chains, variance ratio and efficiency
  validate   end-to-end exactness check: independent kernels must agree
  analysis   analyses of native run directories

Command line: ``fepdmc-toolkit {prepare,run,validate,compare,analyze}``.
Documentation: docs/fepdmc.md.
"""

from __future__ import annotations

FEP_DMC_PIN = "05d08449cffdbd0dfbbbf5009add5cc887bc754b"

from .patches import PatchError  # noqa: E402
from .pooled import side as pooled_side  # noqa: E402
from .prepare import Prepared, check_pin, prepare_tree  # noqa: E402

__all__ = ["FEP_DMC_PIN", "PatchError", "Prepared", "check_pin", "pooled_side", "prepare_tree"]
