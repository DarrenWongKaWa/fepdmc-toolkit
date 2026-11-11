"""Reference IR, interpreter and exact oracles (including the R1 recursive builder) used to
validate the diagram compiler. Not a native drop-in for grouped graphs."""

from .numeric import ABS_ATOL, REL_RTOL, numeric_gate, reference_scale

__all__ = ["ABS_ATOL", "REL_RTOL", "numeric_gate", "reference_scale"]
