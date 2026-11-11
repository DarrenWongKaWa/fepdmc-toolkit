"""R1 grouped-measure toy: a finite signed model with two grouped measures.

  grouped_toy   toy model, grouped measures w_A / w_B, partitions and window operators
  kernels       exact transition kernels and stationarity residuals
  dense_oracle  independent dense enumeration oracle
  run_checks    all checks (``python -m fepdmc_toolkit.r1_grouped.run_checks --out DIR``; needs scipy)

One function of grouped_toy uses fepdmc_toolkit.diagram_compiler (imported lazily).
Derivation: docs/r1-grouped-toy.md.
"""
