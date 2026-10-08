# Changelog

All notable changes to jOSQP are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This file starts with the work that followed upstream commit `36f0613`
("Merge pull request #5"). The first release from this line is **1.0.0**; the
published predecessor is `com.quantego:josqp:0.6.5`.

## [Unreleased]

## [1.0.0] - 2026-10-08

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
  values set to 1 disable these additional guards; they do not imply complete native parity.
  The final inaccurate check retains both guards and the original certificate tolerances.

- **Data and settings ownership**: the public `OSQP` constructor now snapshots both. Reusing
  input `Data` is safe; changing input arrays or settings after construction no longer silently
  mutates an existing solver. Use the explicit update methods. `Model.solve()` now honours
  its current `getParam()` configuration, while `solve(settings)` applies settings each time.
- **Stricter input validation**: malformed CSC arrays, incorrect vector lengths, non-finite
  coefficients, NaN bounds and unsupported runtime options are rejected. One-sided invalid
  bound updates consistently throw `IllegalArgumentException` instead of `IllegalStateException`.
  Empty CSC matrices use an empty index/value array and exactly `n+1` column pointers.
- **Setting types corrected**: `device` is now `int` (default 0) and `cg_tol_reduction` is `int`.
  Nondefault `device`, disabled `allocate_solution` and nonzero `profiler_level` are rejected
  rather than ignored. CG fields remain unused by the direct backend.
- **Automatic iteration cadence**: `adaptive_rho_interval = 0` now selects four times the
  termination interval (100 when termination checking is disabled), matching the iteration
  mode of OSQP 1.0; it is not a timing-based rule.

### Added

- **Index-based bulk API** on `Model.Builder` for large models built in code:
  `addVariables(n, lower, upper)`, `linearObjective(indices, coefficients)`,
  `quadraticObjective(rows, cols, values)` and
  `constraint(indices, coefficients, lower, upper)`, plus
  `Model#setConstraintBounds(int, double, double)`. It avoids boxing, `Pair`
  objects and intermediate maps: a 300,000-variable model with 2,000 rows and
  200 entries per row can be built without per-variable objects or boxed coefficient maps.
  Fresh release measurements: 45.2 ms / 113.4 MB for the 300,000-variable bulk build,
  versus 217.0 ms / 236.1 MB for the algebraic builder, excluding prepared input arrays.
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
  with the supported OSQP 1.0 defaults. This is name coverage, not full C API compatibility:
  Java keeps the boolean `adaptive_rho` iteration mode; time-based and native KKT-based modes
  are not implemented. `polishing` and `warm_starting` are reachable through
  `Settings#isPolishing()/setPolishing(boolean)` and
  `Settings#isWarmStarting()/setWarmStarting(boolean)`; the 0.6 field names
  `polish` and `warm_start` remain for source compatibility.
  `Info` carries `primdual_int` (accumulated `|duality_gap|` at information updates) and
  `rel_kkt_error` (`max(pri_res, dua_res, |duality_gap|)` of the last info
  update), with the Java diagnostic convention documented in README.
- **Benchmark validation**: fresh measurements use the same matrices, settings and limits,
  with one untimed warmup and three independently initialised solver states. The current
  README replaces the earlier mixed-configuration comparisons.
- **Consolidated update entry points** ([issue
  #1](https://github.com/quantego/josqp/issues/1)):
  `OSQP#update_data_vec(q, l, u)` and
  `OSQP#update_data_mat(Px, PxIdx, Ax, AxIdx)` mirror `osqp_update_data_vec` and
  `osqp_update_data_mat`. A `null` argument leaves that vector or matrix
  untouched, one-sided bound updates are allowed, and updating both matrices in
  one call unscales, rescales and refreshes the KKT values once instead of once
  per matrix: on `qplib_8792` (n = 15,129) a combined update takes 8.5 ms against
  17.1 ms for `update_P` followed by `update_A`, the same factor on `qplib_8991`.
  `update_P`, `update_A` and `update_P_A` delegate to it, and
  malformed index arrays (wrong length, out of range, indices without values)
  are rejected instead of silently ignored.
- **Maven release profile** for Maven Central: `mvn -Prelease deploy` attaches
  the sources and javadoc jars, signs all artifacts with GPG and publishes
  through the Sonatype Central Portal. `autoPublish=true`, so a successful
  deployment goes live without a manual release step. `pom.xml` carries the
  metadata Central requires (`name`, `description`, `url`, `licenses`,
  `developers`, `scm`). The published artifact is the plain Java 8 build; the
  Vector API backend stays behind the `vector` profile.
- **GitHub Actions**: CI includes the actual Java 8 test runtime as well as 17/21 and
  optional Vector 26. The release checks tag/POM agreement and serialises publication.
  `release.yml` builds, tests, signs and deploys on every
  `v*` tag; `ci.yml` runs the portable suite on JDK 17 and 21 and the Vector API
  suite on JDK 26 for every push and pull request. The required repository
  secrets and the whole procedure are documented in
  [RELEASING.md](RELEASING.md).

### Changed

- **Test suite migrated to JUnit 5 (Jupiter)**. `mvn test` now drives the suite through the JUnit
  Platform; `junit:junit:4.13.2` is gone and only `hamcrest` remains as an
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

- Preserve unscaled builder data across dimensional and sparsity changes; rebuilt models no
  longer solve a doubly scaled problem. Public solver instances have independent working arrays.
- Compact duplicate CSC triplets after assembly and directly replace existing pattern positions,
  so repeated simulation updates do not retain an unbounded history.
- Sum algebraic cross terms regardless of variable order; bulk quadratic setters retain
  last-assignment semantics. Reject variables from a different model.
- Apply objective offset and direction changes between solves. Offset-only changes need no
  refactorisation; direction changes rebuild P, q and the signed offset together.
- Apply and validate settings on every Model solve, including requested cold starts and
  scaling/sigma changes; combine simultaneous Model P/A updates in one factorisation.
- Preserve infeasibility guards at the iteration limit; insufficient evidence leaves a limit
  status rather than a false inaccurate infeasibility verdict.
- Keep public KKT diagnostics attached to the accepted solution, not a rejected polish candidate.
- Validate combined vector and bulk quadratic updates before any component is changed.
  A failed matrix factorisation marks ERROR and blocks solves until repaired by a valid update.
- Reuse optional linear-residual work arrays across settings updates; disabling the diagnostic
  releases the arrays and stops the checks.


- **The documented terminal command works.** The published jar had no `Main-Class`, so
  `java -jar josqp.jar <file.qps>` failed even though `Parser` has a `main`; the jar plugin now sets
  `Main-Class: com.quantego.josqp.Parser`. Verified with `src/test/resources/qafiro.qps` and against
  the usage message for a missing argument.
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

- Release artifacts now include `LICENSE`, upstream OSQP/QDLDL/AMD `NOTICE`
  and the retained Apache 2.0 / BSD 3-clause license texts in `META-INF`,
  including the sources and Javadoc jars.

- README: the benchmark section now shows nine QPLIB instances measured against
  native OSQP 1.1.3 with the current version (median of three cold starts, both
  configurations spelled out) instead of the historical cross-solver table, which
  only covered an earlier version. The benchmark harness and the native comparison
  tooling live outside the repository.
- README: quadratic-term convention and the migration note for existing models,
  how model updates reuse or rebuild the solver, the settings that OSQP 1.0 adds
  and the ones jOSQP cannot honour, the index-based bulk API, the infeasibility
  certificates section, building from source, the license and a status note on
  the addressed issues. The code examples were compiled against the released jar
  with `--release 8`.
- `JOSQP-Vector-API.md` (project report) documents the measurements, the native
  OSQP/PIQP comparisons, the GitHub issue assessment and the remaining OSQP 1.0
  API delta.
