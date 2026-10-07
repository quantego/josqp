# Changelog

All notable changes to jOSQP are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This file starts with the work that followed upstream commit `36f0613`
("Merge pull request #5"). The first release from this line is **1.0.0**; the
published predecessor is `com.quantego:josqp:0.6.5`.

## [Unreleased]

## [1.0.0] - 2026-10-06

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
  objects and intermediate maps: a 300,000-variable model with 2,000 rows and
  200 entries per row builds in about 41 ms and allocates 129 MB, where the
  algebraic builder path needs 224 ms and 236 MB. That path itself came down
  from 322 ms, 343 MB and 694 MB peak heap to 224 ms, 236 MB and 524 MB peak
  through the primitive triplet storage.
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
- **OSQP 1.0 settings parity** ([issue #1](https://github.com/quantego/josqp/issues/1)):
  `device`, `linsys_solver`, `allocate_solution`, `profiler_level`, `rho_is_vec`,
  `cg_max_iter`, `cg_tol_reduction`, `cg_tol_fraction` and `cg_precond` now exist
  with the OSQP 1.0 defaults, so code written against `OSQPSettings` compiles
  unchanged. `polishing` and `warm_starting` are reachable through
  `Settings#isPolishing()/setPolishing(boolean)` and
  `Settings#isWarmStarting()/setWarmStarting(boolean)`; the 0.6 field names
  `polish` and `warm_start` remain for source compatibility.
  `Info` carries `primdual_int` (running integral of `|duality_gap|`) and
  `rel_kkt_error` (`max(pri_res, dua_res, |duality_gap|)` of the last info
  update), both filled like in OSQP 1.0.
- **Maven release profile** for Maven Central: `mvn -Prelease deploy` attaches
  the sources and javadoc jars, signs all artifacts with GPG and publishes
  through the Sonatype Central Portal. `autoPublish=true`, so a successful
  deployment goes live without a manual release step. `pom.xml` carries the
  metadata Central requires (`name`, `description`, `url`, `licenses`,
  `developers`, `scm`). The published artifact is the plain Java 8 build; the
  Vector API backend stays behind the `vector` profile.
- **GitHub Actions**: `release.yml` builds, tests, signs and deploys on every
  `v*` tag; `ci.yml` runs the portable suite on JDK 17 and 21 and the Vector API
  suite on JDK 26 for every push and pull request. The required repository
  secrets and the whole procedure are documented in
  [RELEASING.md](RELEASING.md).

### Changed

- **Test suite migrated to JUnit 5 (Jupiter)**. `mvn test` and `benchmarks/run.py` now drive the
  suite through the JUnit Platform; `junit:junit:4.13.2` is gone and only `hamcrest` remains as an
  extra test dependency for the two matcher assertions. JUnit 5 expects the failure message as the
  *last* argument, so every message-carrying assertion was reordered. The archived `baseline` and
  `pre-workspace` harness modes keep running the historical JUnit 4 sources.
- Settings that jOSQP cannot honour are now rejected instead of silently
  ignored: `linsys_solver = OSQP_INDIRECT_SOLVER` (no conjugate-gradient solver)
  and `rho_is_vec = false` (scalar rho) throw an `IllegalArgumentException` from
  the `OSQP` constructor. `device`, `allocate_solution`, `profiler_level` and
  the four CG parameters are inert for a pure Java solver; the
  `SettingsConformanceTest` documents them as such.
- `OSQP` validates its settings in the constructor, as `osqp_setup` does in C,
  so an invalid configuration fails at construction instead of at the first
  solve.
- `QDLDL` keeps values in locals instead of reading and writing the work array
  repeatedly, and the diagonal scaling is folded into the backward solve.
- `CSCMatrixBuilder` collects primitive `(row, column, value)` triplets and
  assembles CSC with a stable counting sort instead of one `TreeMap` per column;
  `Model` stores `q`, `l`, `u` and `rhs` in primitive lists.
- `OSQP.Data#takeOwnership` hands freshly built arrays to the solver without a
  second copy.

### Fixed

- **`mvn test` compiles again.** The pom depended on
  `org.hamcrest:hamcrest-library:2.2`, which is only a shim without
  `org.hamcrest.Matchers`; the tests compiled through the project's own
  `benchmarks/run.py` but never through Maven. The dependency is now
  `org.hamcrest:hamcrest:2.2` with JUnit's transitive `hamcrest-core:1.3`
  excluded. Verified with Maven 3.10.0 on JDK 26: 90 tests, plus 93 with
  `-Pvector`.
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

### Documentation

- README: quadratic-term convention and the migration note for existing models,
  how model updates reuse or rebuild the solver, and the infeasibility
  certificates section.
- `JOSQP-Vector-API.md` (project report) documents the measurements, the native
  OSQP/PIQP comparisons, the GitHub issue assessment and the remaining OSQP 1.0
  API delta.
