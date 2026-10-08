# jOSQP

jOSQP is a fork of the quadratic programming solver [OSQP](http://osqp.org) (Operator Splitting
Quadratic Program), written entirely in Java. It has no runtime dependencies beyond Java 8, and its
supported settings and termination criteria follow OSQP 1.0, with the differences described below.

The solver can handle quadratic programming problems of the form: 

$$
\begin{align}
&& \min \ & \frac{1}{2}x^\top Px + q^\top x + c \\
&& \mathrm{s.t.} \ & l \leq Ax \leq u
\end{align}
$$

where $x\in\mathbb{R}^n$ represents variables, $A$ is an $m \times n$ matrix of constraints with lower bounds $l$ and upper bounds $u$, and $P$ is a positive semi-definite matrix ($P\in S^n_+$).

The solution algorithm of the OSQP solver is based on the ADMM (Alternating Direction Method of Multipliers) which is described [in this paper](https://arxiv.org/abs/1711.08013).

## Usage

While OSQP has been linked with a few Python packages to facilitate model creation, jOSQP provides its own algebraic model builder as well as a parser for MPS and QPS files.

### Model Builder

The jOSQP builder class allows creating variables, constraints, and linear expressions as Java objects. Unfortunately, there is no operator overloading in Java, so that use of mathematical operators is not possible. Nonetheless, decision variables as well as their bounds, constraint, and objective function can be created in a single line. For example, a simple quadratic program such as

$$
\begin{align}
&&\min\ & 5x + 0.5x^2 + 3y \\
&&\mathrm{s.t.}\  & 2x-y \geq 4 \\
&& &1 \leq  x \leq 3 \\
&& &0 \leq y \leq 5
\end{align}
$$

can be expressed and solved in only a few lines of Java code.

```
import com.quantego.josqp.Model;
import com.quantego.josqp.Model.Variable;

Model.Builder builder = Model.getBuilder();
Variable x = builder.addVariable().lb(1).ub(3);
Variable y = builder.addVariable().ub(5);
builder.setObjective().add(5, x).add(0.5, x, x).add(3, y).minimize();
builder.addConstraint().add(2, x).add(-1, y).geq(4);
Model model = builder.build();
model.solve();
```

`Variable`, `Builder`, `Objective` and `Constraint` are nested types of `Model`; the model is built
and solved with Java 8, no language level above that is required.

### Bulk model API

The algebraic builder creates one `Variable` object and one map entry per coefficient, which becomes
expensive for models with hundreds of thousands of entries. `Builder` therefore also offers an
index-based API that works on primitive arrays, and `Model#setConstraintBounds` changes bounds later
without touching the sparsity pattern:

```
int n = 300_000;
double[] lower = new double[n];
Model.Builder builder = Model.getBuilder();
builder.addVariables(n, lower, null); // null means no upper bound
builder.linearObjective(new int[] {0, 1}, new double[] {1.0, 2.0});
builder.quadraticObjective(new int[] {0}, new int[] {0}, new double[] {0.5});
builder.constraint(new int[] {0, 1}, new double[] {2.0, -1.0}, 4.0, Double.POSITIVE_INFINITY);
Model model = builder.build();
```

The bulk path avoids per-variable objects and coefficient maps. Matrix value updates replace
existing entries, so a simulation with a fixed sparsity pattern does not retain its write history.

### Quadratic terms

`add(a, x, x)` sets the coefficient of `x*x`, `add(a, x, y)` the coefficient of `x*y`, and
`quadraticObjective(rows, cols, values)` does the same for index pairs. The coefficient is the
coefficient of the monomial, unchanged by the argument order, and the entry is always stored in the
upper triangle of `P`.

Earlier versions doubled the off-diagonal coefficient, so `add(-1, x, y)` used to mean `-2*x*y` and
`add(-1, y, x)` raised `quad_form matrix is not upper triangular` ([issue
#2](https://github.com/quantego/josqp/issues/2)). Models that relied on the halving workaround
(`add(-0.5, x, y)` for `-x*y`) have to pass the full coefficient now. This is a breaking change;
[CHANGELOG.md](CHANGELOG.md) lists it together with the other changes since upstream `36f0613`.

### Updating a model

A `Model` can be solved repeatedly. Changing an objective or constraint coefficient between solves
reuses the current factorisation. Adding or removing a nonzero changes the sparsity pattern; jOSQP
then rebuilds the solver and runs a new symbolic factorisation, because the KKT mappings and the
scaling are tied to the pattern. At the solver level `OSQP.update_A` and `OSQP.update_P` stay
value-only, exactly like `osqp_update_data_mat` in OSQP: they require the pattern to be unchanged
([issue #4](https://github.com/quantego/josqp/issues/4)).

### Settings

`new OSQP(data, settings)` copies and validates the unscaled problem data and settings before
creating its workspace. A `Data` instance can be reused by several independent solvers. Mutating
that instance afterwards does not change an existing solver; use the update methods instead.
`Data` itself copies the arrays supplied to its public constructor. `Model` uses an internal
ownership-transfer path for freshly allocated arrays while preserving its unscaled builder data.

`Model.solve(settings)` applies the supplied settings on every call. Runtime changes retain the
workspace; changes to scaling, initial cost scaling or sigma rebuild it. `Model.solve()` uses the model's current
`getParam()` configuration, and `solve(null)` requests fresh defaults. The direct solver's
`update_settings(settings)` accepts runtime changes, but rejects scaling/sigma/initial-cost changes which
require a new solver.

The supported numerical defaults match OSQP 1.0: `eps_abs = eps_rel = 1e-3`, `max_iter = 4000`,
`rho = 0.1`, `adaptive_rho_interval = 50`, `check_dualgap = true`, `time_limit = 1e10`, polishing
off and warm starting on. This is a pure Java implementation of the direct QDLDL backend, not a
complete reproduction of the C runtime API:

- `polish` and `warm_start` retain their older Java field names. OSQP 1.0 names are available through
  `setPolishing` / `isPolishing` and `setWarmStarting` / `isWarmStarting`.
- `adaptive_rho` retains a boolean: enabled means iteration-based updates. Interval `0` selects
  four times `check_termination`, or 100 when termination checking is disabled. The C API's time-
  and KKT-based mode enum is not implemented. `adaptive_rho_kkt` is an experimental Java progress gate.
- `device` and `cg_tol_reduction` are integers like the OSQP 1.0 fields. Only `device = 0`,
  `allocate_solution = true` and `profiler_level = 0` are supported; other values are rejected.
- The indirect solver and `rho_is_vec = false` are rejected. CG settings are validated but unused
  with the direct solver. `OSQP_UNKNOWN_SOLVER` resolves to the direct backend.
- `Info.primdual_int` accumulates absolute gap values at information updates, rather than elapsed
  time. `Info.rel_kkt_error` is `max(pri_res, dua_res, abs(duality_gap))` for the accepted solution.
  Rejected polishing candidates do not replace these solution diagnostics.

The low-level `update_data_vec(q,l,u)` and `update_data_mat(Px,PxIdx,Ax,AxIdx)` methods allow `null`
for unchanged data. Matrix indices are positions in the existing CSC value array. Invalid lengths,
indices, NaN values and contradictory bounds are rejected before replacement. A numerical matrix
factorisation failure marks `ERROR` and blocks subsequent solves until a successful valid matrix
update; it does not roll back the attempted matrix update. `Model` rebuilds automatically when the
sparsity pattern changes, and submits simultaneous P/A value changes in one combined update.

### Infeasibility certificates

Infeasibility detection uses the most recent primal/dual step as a numerical certificate. Loose
certificate tolerances can mistake a nearly null direction for evidence of infeasibility, even on
feasible poorly scaled problems. Two Java defaults make this test more conservative:

- `infeasibility_streak = 2`: require consecutive checks on distinct iterates.
- `infeasibility_margin = 0.1`: multiply the requested certificate tolerance by this safety factor.

Setting both to `1` disables the additional guards; it does not guarantee identical behaviour to
native OSQP in every case. The final approximate check relaxes optimality tolerances by ten, but
never relaxes certificate tolerances or bypasses persistence. A weak certificate at the iteration
limit therefore leaves `MAX_ITER_REACHED`. The guards reduce false verdicts; they are heuristics,
not a proof of feasibility or infeasibility. Check residuals and status when accuracy matters.

The optional `check_linear_residual`, `linear_refine_iter` and `adaptive_rho_kkt` diagnostics are
experimental and disabled by default. Nonconvex QPs are outside the solver's supported problem class.

### Difficult convex problems

A limit status does not establish infeasibility. Automatic residual balancing is a heuristic;
for some badly conditioned problems it chooses rho values which slow or destabilise convergence.
Keep the gap check enabled and try a fixed rho with a sufficient iteration budget:

```java
OSQP.Settings settings = new OSQP.Settings();
settings.adaptive_rho = false;
settings.rho = 0.1;                 // use 10.0 for Maros-Meszaros primalc5
settings.max_iter = 100000;
settings.time_limit = 120;
settings.eps_abs = settings.eps_rel = 1e-9;
settings.eps_prim_inf = settings.eps_dual_inf = 1e-12;
settings.polish = true;
model.solve(settings);
```

On the unchanged Maros-Meszaros input matrices, this profile solves `primalc1`, `primalc2` and
`primalc8` at rho 0.1 and `primalc5` at rho 10. Their original-matrix feasibility, stationarity
and gap checks pass; native OSQP and PIQP independently agree on their objectives. Polishing
is accepted only when it preserves the requested solution tolerances, including the gap.
These are tested problem-specific settings, not new defaults or a general convergence guarantee.
Increasing a time limit alone is not a substitute for examining the residuals.

Two additional Java conditioning settings leave OSQP's defaults intact:
`scaling_cost = 1` is an initial positive objective multiplier before Ruiz equilibration;
`rho_eq_over_rho_ineq = 1000` weights equalities relative to inequalities. Costs, duals and
residuals are reported in the original units. Changing `scaling_cost` requires rebuilding;
`Model` handles that, and ratio changes refactor the existing workspace. Non-unit cost scaling
requires `scaling > 0`; unrepresentable scaling and overflowing rho ratios are rejected.

QPLIB 10038 (160,800 variables) is solved and independently validated with this advanced profile:

```java
settings.eps_abs = settings.eps_rel = 1e-6;
settings.rho = 0.1;
settings.scaling_cost = 1000;
settings.rho_eq_over_rho_ineq = 1e6;
settings.max_iter = 200000;
settings.delta = 1e-8;
settings.polish_refine_iter = 20;
model.solve(settings);  // adaptive_rho=false, polish=true, time_limit=120 as above
```

Its objective agrees with PIQP at approximately -0.05167742451254; the measured maximum bound
violation is 1.4e-13. These results do not establish a universal parameter choice. QPLIB 10034
remains unresolved at strict accuracy in the diagnostic runs. At a time limit, JOSQP now checks
both exact and ten-times relaxed optimality tolerances before reporting `TIME_LIMIT_REACHED`;
qualified iterates are explicitly marked `SOLVED_INACCURATE` when needed.

### MPS/QPS Reader

Also, jOSQP comes with a parser for MPS and QPS files to import and solve existing linear and quadratic programming problems. Use the static function `Parser.readQmps` to read a file and pass the resulting model data to the solver:
```
Parser p = Parser.readQmps(<qps_file_name>);
OSQP opt = new OSQP(p.getData(), new OSQP.Settings());
opt.solve();
```
The jar is executable and solves a file from the terminal:
```
java -jar josqp.jar <qps_file_name>
java -jar josqp.jar src/test/resources/qafiro.qps
```
After solving the problem through the terminal, jOSQP creates a log file, `josqp.log`, that contains a summary of the solution process. The terminal command has been used for benchmarking.

## Benchmarks

### Native OSQP comparison

Remeasured on 2026-10-08 after the short-column optimization, using portable solver commit
`62d29b4`: nine QPLIB problems plus Maros-Meszaros `liswet3` and `ubh1`, on identical exported,
unscaled matrices, macOS arm64 / OpenJDK 26.0.2. Both solvers use `eps_abs = eps_rel = 1e-6`,
`eps_prim_inf = eps_dual_inf = 1e-8`, iteration-based adaptive rho with interval 50, 100,000
iterations, a 120-second limit, gap checking and no polishing. The extra Java infeasibility
guards are disabled for this comparison. The Java artifact uses no Vector API.

Each solver gets one untimed warmup and three measured solves with independently initialised
solver states; these are **warm JVM measurements, not three cold JVM processes**. Solver runs
are sequential. The table shows median solve wall time, excluding parsing, export and setup.
Native OSQP 1.1.3 runs its C core through the Python wrapper; the measured call includes that
wrapper's overhead.

| Instance | n | Java | native OSQP | Java / native | Status |
|---|---:|---:|---:|---:|---|
| `qplib_8845` | 1,546 | 7.029 s | 7.839 s | 0.897 | SOLVED / solved |
| `qplib_9002` | 2,890 | 9.142 s | 9.732 s | 0.939 | SOLVED / solved |
| `qplib_8938` | 4,001 | 12.922 s | 14.519 s | 0.890 | MAX_ITER_REACHED / maximum iterations reached (excluded) |
| `qplib_8906` | 5,223 | 0.945 s | 1.033 s | 0.915 | SOLVED / solved |
| `qplib_8559` | 10,000 | 31.749 s | 35.132 s | 0.904 | SOLVED / solved |
| `qplib_8616` | 13,870 | 0.056 s | 0.052 s | 1.071 | SOLVED / solved |
| `qplib_8991` | 14,400 | 0.193 s | 0.184 s | 1.044 | SOLVED / solved |
| `qplib_8792` | 15,129 | 0.378 s | 0.369 s | 1.024 | SOLVED / solved |
| `qplib_8602` | 34,552 | 32.196 s | 32.008 s | 1.006 | SOLVED / solved |
| `liswet3` | 10,002 | 10.911 s | 13.838 s | 0.788 | SOLVED / solved |
| `ubh1` | 18,009 | 8.085 s | 8.644 s | 0.935 | SOLVED / solved |

All ten jointly solved cases pass independent primal, stationarity and gap checks on the original
data and agree in objective within relative tolerance 1e-5 (with an absolute floor of 1e-5).
`qplib_8938` reaches 100,000 iterations on both sides and is excluded from the geometric means.
The ratios below are Java / native; a value below 1 means less time for Java.

| Validated subset | Cases | Geometric mean solve | Geometric mean setup + solve |
|---|---:|---:|---:|
| QPLIB | 8 | 0.973 | 1.019 |
| QPLIB + the two Maros-Meszaros cases | 10 | 0.949 | 0.985 |

The expanded sample uses about 5% less solve time for Java and is approximately equal including
setup. `liswet3` takes about 21% less solve time; `ubh1` about 6% less. The eight-case QPLIB sample
remains approximately equal. This is a selected set, not a full-suite comparison or a general
Java speed advantage. Small differences depend on runtime and measurement variation.

Setup includes solver validation/copying, scaling and initial factorisation; parsing, export and
preparation of caller-owned arrays or SciPy matrices are excluded. Setup + solve uses the median
of each run's combined time, rather than adding separate medians. These fresh figures supersede
the preceding Java/native table; the Java A/B measurements below use their separately stated profile.

Full runs, setup times, iterations, original-matrix quality checks, input checksums and build
provenance are retained in `release-current-bench.json` and `release-current-summary.json` in the
external `josqp-tools` measurement workspace. No Python code or runtime dependency is required
by this Java repository.

### Portable short-column optimization

The direct solver selects an unrolled triangular-solve kernel when every factor column has
at most four off-diagonal entries. It makes the short loops explicit for the JIT, preserves
the original arithmetic order and requires no Vector API or newer Java language features.
Wider factors keep the general kernel. Selection happens once per factor, without a plan array
or per-iteration allocations.

Matched Java A/B measurements on macOS arm64, OpenJDK 26.0.2, against the immediately preceding
release candidate (`3b4f5fe`), excluding parsing and setup:

| Instance | Previous solve | Short-column selection | Relative time saved |
|---|---:|---:|---:|
| `liswet3` | 15.974 s | 11.047 s | 30.8% |
| `ubh1` | 9.830 s | 8.181 s | 16.8% |
| `qplib_8991` | 0.190 s | 0.190 s | approximately unchanged |
| `primalc8` | 0.427 s | 0.427 s | approximately unchanged |
| `qplib_8559` | 32.190 s | 32.120 s | approximately unchanged |
| `qplib_8602` | 32.312 s | 32.276 s | approximately unchanged |

All measured solves return `SOLVED`, with bitidentical full primal/dual vectors and iteration
counts between variants; independent original-matrix feasibility, stationarity and gap checks
pass. The measured heap allocations during setup and solve are unchanged. Code cache, metaspace
and peak RSS were not measured. The geometric mean solve-time ratio over these **six selected
instances** is 0.912 (about 9% less time); this is not a full-suite or native-solver comparison.

`liswet3` and `ubh1` each use three new JVMs per variant, alternating A/B order, with one untimed
warmup and two measured fresh solver states in each JVM. The four other cases each use one JVM
per variant with one warmup and three measured states. Both variants use eps_abs=eps_rel=1e-6,
eps_prim_inf=eps_dual_inf=1e-12, default certificate guards, rho interval 50, max_iter=1,000,000,
time_limit=120 and no warm starting or polishing. `primalc8` instead uses eps_abs=eps_rel=1e-9,
fixed rho=0.1 and polishing. These are warm-JIT measurements; gains depend on the factor pattern
and runtime. ADMM loop fusion and direct RHS/permutation fusion were tested and not retained
because they did not show consistent improvements.

### Model construction

Remeasured on 2026-10-08 with the current portable solver: median of three new JVM processes
per size and mode, OpenJDK 26.0.2, `-Xmx2g`; mode order alternates between repetitions.
The same primitive input arrays are prepared before timing either builder; measured allocation
excludes those arrays. Construction
ends at `Model.build()`, before solver setup, CSC assembly and factorisation.

| Variables | User rows × entries | Algebraic build / allocation | Bulk build / allocation |
|---:|---:|---:|---:|
| 100,000 | 1,000 × 200 | 85.6 ms / 85.1 MB | 21.0 ms / 36.8 MB |
| 300,000 | 2,000 × 200 | 210.4 ms / 236.1 MB | 46.9 ms / 113.4 MB |

Here MB means allocated bytes / 1,000,000. A fixed-pattern stress test performs 100,000 rounds of
three coefficient assignments and retains exactly three current entries. Solver rebuilds retain
unscaled builder data and allocate independent working matrices; this ownership is essential for
correct solutions and adds storage compared with sharing those arrays unsafely.

The optional Vector backend remains experimental. Sparse factor workspaces and combined matrix
updates are retained; the default portable path uses neither native libraries nor incubator APIs.

## Installation

The planned 1.0.0 Maven Central coordinates are shown below. Until that release is published,
build locally or use the previously published version. After publication, add this dependency:
```
<dependency>
  <groupId>com.quantego</groupId>
  <artifactId>josqp</artifactId>
  <version>1.0.0</version>
</dependency>
```

Another option is to simply download the [latest build](https://github.com/quantego/josqp/releases/latest) of the jar from the release page and import it in your Java project.

## Building from source

The build uses Maven and produces a Java 8 artifact (`maven.compiler.release=8`):

```
mvn test                            # portable suite, JUnit 5
mvn -Pvector test                   # needs JDK 26 and jdk.incubator.vector
mvn -Prelease package               # main, sources and javadoc jar in target/
java -jar target/josqp-1.0.0.jar src/test/resources/qafiro.qps
```

## Requirements

* Java 8 or higher to use the library; the published jar is compiled with `--release 8`
* JDK 26 for the optional Vector API profile
* A JDK and Maven 3.9 or newer to build and test

## Status

Version 1.0.0 fixes [issue #2](https://github.com/quantego/josqp/issues/2) (quadratic coefficients)
and [issue #4](https://github.com/quantego/josqp/issues/4) (model updates), and includes the supported
OSQP 1.0 settings and consolidated vector/matrix update calls from
[issue #1](https://github.com/quantego/josqp/issues/1). The remaining C API differences are listed
above. [CHANGELOG.md](CHANGELOG.md) includes migration notes and breaking changes.

## License

The jOSQP-specific work uses the MIT License, see [LICENSE](LICENSE).
The port retains upstream OSQP/QDLDL (Apache 2.0) and AMD (BSD 3-clause) terms;
see [NOTICE](NOTICE) and [THIRD-PARTY-LICENSES](THIRD-PARTY-LICENSES).
All three release jars include these texts under `META-INF`.

## Authors
* Nils Loehndorf
* Farid Alavi
