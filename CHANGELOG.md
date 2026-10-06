# Changelog

All notable changes to jOSQP are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This file starts with the work that followed upstream commit `36f0613`
("Merge pull request #5"). Everything below is **not yet released**; `pom.xml`
still declares `0.6.6`.

## [Unreleased]

### Changed — breaking

- **Quadratic objective terms now use the coefficient of the monomial**
  ([issue #2](https://github.com/quantego/josqp/issues/2)).
  `Objective.add(a, x1, x2)` and `Builder.quadraticObjective(...)` set the
  coefficient of `a * x1 * x2`. Up to now the off-diagonal contribution was
  doubled, so `add(-1, x, y)` built the term `-2*x*y` and the documented `-x*y`
  required the workaround `add(-0.5, x, y)`.
  - `add(a, x, x)` is unchanged: it stores `a` as the coefficient of `x*x`.
  - Code that used the halving workaround must pass the full coefficient now,
    otherwise the cross term is halved and the optimum moves.
  - Code that wrote the coefficient without the workaround silently optimised a
    different, usually unbounded problem before; for
    `min x^2 + y^2 - x*y - x - y` jOSQP reported `DUAL_INFEASIBLE` and now
    returns `SOLVED`, (1, 1), objective -1.
  - The entry is always stored in the upper triangle of `P`, independent of the
    argument order: `add(a, y, x)` no longer raises
    `quad_form matrix is not upper triangular`.
  - `Model.setQbj(int, int, double)` follows the same convention.

- **Default settings aligned with OSQP 1.0**: `adaptive_rho_interval` 0 -> 50,
  `RHO_EQ_OVER_RHO_INEQ` 1e4 -> 1e3, disabled `time_limit` +inf -> 1e10.
  These change the ADMM path and therefore iteration counts and timings.

- **Infeasibility verdicts are more conservative**: a certificate must now hold
  over two consecutive termination checks
  (`Settings#infeasibility_streak`, default 2) and has to clear the tolerance by
  an order of magnitude (`Settings#infeasibility_margin`, default 0.1). Both
  values set to 1 reproduce the upstream OSQP verdicts exactly.

### Added

- **Index-based bulk API** on `Model.Builder` for large models built in code:
  `addVariables(n, lower, upper)`, `linearObjective(indices, coefficients)`,
  `quadraticObjective(rows, cols, values)` and
  `constraint(indices, coefficients, lower, upper)`, plus
  `Model#setConstraintBounds(int, double, double)`. It avoids boxing, `Pair`
  objects and intermediate maps; building a 300,000-variable model takes about
  41 ms instead of 322 ms and allocates 129 MB instead of 343 MB.
- **Optional Vector API backend** (`src/vector`, Maven profile `vector`,
  `-Djosqp.vector=true`). It accelerates isolated norm computations by up to
  7.5x, but the contribution to a full solve stays below 1 %, so the default
  artifact remains plain Java 8 compatible. The portable test suite passes on
  JDK 25 and 26.
- `OSQP.Settings#check_dualgap`, the OSQP 1.0 duality-gap termination criterion
  (enabled by default like in OSQP 1.0).
- `OSQP.Settings#infeasibility_streak` and
  `OSQP.Settings#infeasibility_margin` (see above).
- `Model#getNumVariables()` and `Model#getNumConstraints()`.

### Fixed

- **Repeated solves after changing matrix values**
  ([issue #4](https://github.com/quantego/josqp/issues/4)).
  `Model#solve` passed CSC *row indices* to `OSQP#update_A`/`update_P`, which
  expect element *positions*. Unless the sparsity pattern happened to run along
  the diagonal - which the structural variable-bound rows in `A` and any
  off-diagonal entry in `P` prevent - the values landed in the wrong entries.
  The result was a silently wrong optimum, or
  `Invalid KKT inertia: expected 2 positive pivots, got 1` for a changed cross
  term. A value update is applied in place only while the pattern is unchanged;
  adding or removing a nonzero now rebuilds the solver with a new symbolic
  factorisation and warm-starts it from the previous solution.
- **Feasible problems no longer report infeasibility at the default tolerances.**
  `hues-mod` is solved in 875 iterations (objective 3.4761704979875020e7) and
  `huestis` in 6675 iterations (3.4755179890390e11), matching native OSQP at
  `eps_prim_inf = eps_dual_inf = 1e-5`. `primalc5` now ends in
  `MAX_ITER_REACHED` instead of a wrong `DUAL_INFEASIBLE`. Genuine
  infeasibility problems are still detected at the same iteration.
- **Parser**: BOUNDS entries without a value (`PL bounds x2`) no longer read the
  bound name as a variable name; `RANGES` sections are parsed and validated
  correctly; MPS names are treated as literals instead of regular expressions;
  invalid `QUADOBJ` entries are rejected.
- **CSC assembly, polishing and KKT pivot check**: sorting, the polishing
  acceptance test and the negative-pivot detection were corrected.
- **Resource handling**: parser exceptions are propagated instead of swallowed,
  and the file reader is closed on parse errors.
- **Time limit**: it is checked monotonically and independent of the logging and
  termination-check interval; `0` disables it as documented.
- **Updates**: repeated solves no longer inherit stale status or timing;
  bound updates handle a switch between equality and inequality rows in the rho
  update; an invalid single-bound update no longer leaves the solver partly
  modified; unbounded constraints no longer produce `0 * Infinity` NaNs in the
  infeasibility test.
- **Factorisation**: the symbolic LDL analysis, the elimination schedule and the
  factor/work arrays are reused across rho updates. For `qplib_8602` this
  reduces the solve time by 13.9 % (31.56 s -> 27.17 s) and the temporary
  allocations from 767 MB to 64 bytes.

### Changed

- `QDLDL` keeps values in locals instead of reading and writing the work array
  repeatedly, and the diagonal scaling is folded into the backward solve.
- `CSCMatrixBuilder` collects primitive `(row, column, value)` triplets and
  assembles CSC with a stable counting sort instead of one `TreeMap` per column;
  `Model` stores `q`, `l`, `u` and `rhs` in primitive lists.
- `OSQP.Data#takeOwnership` hands freshly built arrays to the solver without a
  second copy.

### Documentation

- README: quadratic-term convention and the migration note for existing models,
  how model updates reuse or rebuild the solver, and the infeasibility
  certificates section.
- `JOSQP-Vector-API.md` (project report) documents the measurements, the native
  OSQP/PIQP comparisons, the GitHub issue assessment and the remaining OSQP 1.0
  API delta.
