# Open directions

Not implemented. Listed so that a contributor knows where the code stops.

- **Grouped Monte Carlo on a real material.** The R1 grouped toy
  ([r1-grouped-toy.md](r1-grouped-toy.md)) shows that a grouped measure
  built from summed window operators can be sampled exactly, and the diagram
  compiler can supply the window operator. On the toy the sign gain is below
  10%, and the fixed-order R1 evaluator showed no speed advantage over CSE
  ([results.md](results.md)). A material version would need external phonon
  legs, larger windows, a native implementation and an efficiency test
  before it is worth building.
- **Variable-order compilation.** The compiler handles fixed-order families
  only; a variable-order sampler needs a measure, proposal and reverse maps
  and a normalization contract ([diagram-compiler.md](diagram-compiler.md)).
- **Learned stage-1 scores.** P1's cheap score is analytic. A learned
  predictor of the residual `ell_R − ell_prop` could sharpen stage 1 while
  stage 2 keeps the target exact. No data were collected and no model was
  trained.
- **P1 beyond the current domain.** Finite temperature, multiband,
  `sample_gt` and other DMC methods are refused today; each needs a
  target-equivalence fixture before it can be enabled.
- **Converged formation energies.** The LiF numbers are fixed-setting
  benchmarks. A 1%-level ground-state value would need the projection,
  order and numerical remainders closed and many more chains.
