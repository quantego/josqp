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

Full companion-suite census measured on 2026-10-08 using portable solver commit
`c455fde` (numerical classes unchanged from `62d29b4`), macOS arm64 /
OpenJDK 26.0.2, against the builtin C solver of OSQP Python package 1.1.3. The companion QPS inputs cover all 19
continuous convex QPs with linear or box constraints in the
[QPLIB catalogue](https://qplib.zib.de/instances.html), plus 139 available
Maros-Meszaros files/variants. QSCORPIO was recovered from the
[Maros-Meszaros mirror](https://github.com/optimizers/maros-meszaros-mirror) after discovering
a broken companion archive. The listed `xdata.qps` input is missing; it is not classified as
infeasible. No synthetic collection was present in the companion repository.

Each pair uses the same exported original matrices and numerical settings: `eps_abs =
eps_rel = 1e-6`, certificate tolerances `1e-12`, scaling 10, rho 0.1, gap checking,
no polishing, and zero initial iterates in fresh solver states. A means iteration-based adaptive rho with interval 50;
F means fixed rho. The four `primalc` problems use the previously validated fixed-rho
profile with `eps_abs = eps_rel = 1e-9`; each such profile is applied to both solvers.
Failed screenings additionally try the same fixed-rho profiles, 0.1 then 10, on
both sides (60-second probes, extended to 225 seconds if one solver validates a solution);
validated profiles are then measured with the full 225-second budget. `qshare1b` uses `1e-7`, and QPLIB 10034/10038 use `1e-9` to test their small
objectives. The first jointly validated profile is measured, rather than choosing a
separate fastest configuration for each solver. Profiles are identified in the tables;
this compares configured solvers, not unchanged default settings.
Java's additional certificate guards are disabled. The portable artifact uses no Vector API.

The initial budget is 180 seconds and 1,000,000 iterations per solve; unsuccessful
limit-bound pairs receive a **25% extension to 225 seconds and 1,250,000 iterations**.
Workers are sequential, with a 4 GiB Java heap and a 4.5 GiB process-RSS cap on both sides.
Resource failures and exhausted budgets remain visible below and are not called infeasible.
The process cap also includes parsing/export and adapter memory, so a resource failure does
not establish whether the numerical core could solve the instance on a larger machine.

Measured pairs get at least one second of untimed warmup, followed by independently
initialised solver states: 100 repetitions for sub-5-ms screenings, 15 below 20 ms and
three otherwise. Each measured pair uses one JVM/native-worker process; this sweep does
not estimate variation across fresh JVM forks or alternate solver order. These are warm-JVM timings. Java sets `warm_start=false`; native retains
its `warm_starting=true` default, but every call receives a newly created solver, never a
previous solution. Native clears initial iterates in setup; Java explicitly clears them
at solve entry. This boundary difference matters most for tiny problems; the combined
setup + solve column includes both phases. Solve wall time excludes parsing, export and
setup; native timings include the Python wrapper call around the compiled C solver.
Setup includes validation/copying, scaling and initial factorisation; caller-side SciPy
matrix construction is excluded. Setup + solve is the median of each run's combined time.

**The raw setup + solve ratio is not matched work.** The native adapter takes its setup timestamp
in Python around `solver.setup()` and its solve timestamp around `solver.solve()`, while the Java
runner uses `System.nanoTime()` around `new OSQP(...)` and `OSQP.solve()`. Marshalling the SciPy
matrices into the C solver happens *inside* `solver.setup()`, so it is charged to native setup even
though the Java side has no counterpart: measured against OSQP's own `info.setup_time`, which starts
inside the compiled solver, that adapter share is a median of
32&nbsp;% of the recorded native setup time and reaches 53&nbsp;% on the largest
matrices. The `setup + solve (adapter overhead removed)` column subtracts the measured adapter share
per pair; the measurement script and the raw data live outside this Java repository in
`josqp-tools/review/setup_wrapper_overhead.py` and `josqp-tools/results/setup-wrapper-overhead.json`.
Solve times need no correction: over all pairs
below 50 ms native solve time the measured adapter share of `solver.solve()` is not positive
(median -0.12 ms), and above that it is below one percent. The geometric mean of the solve column is
therefore already matched work.

A pair enters the geometric means only when every measured solution passes independent
original-matrix feasibility, stationarity and gap checks at the requested tolerance, and
objectives agree within `10 * eps * max(1, abs(objective))`. `SOLVED_INACCURATE` alone is
insufficient. QPLIB 10038 additionally must agree with its independent PIQP reference;
10034 must improve on the previously recorded approximate PIQP candidate (which is not
claimed as a proven optimum). This prevents weak absolute tolerances from admitting
near-zero outcomes inconsistent with the independent diagnostic runs.
Ratios are Java / native; values below one favour Java. No small-instance cutoff is used. Identical matrix exports (including `hs268`/`s268`)
are listed under both input names but counted only once in the geometric means.

| Collection | Available inputs | Validated distinct pairs | Geometric mean solve | Geometric mean setup + solve (measured) | setup + solve (adapter work removed) |
|---|---:|---:|---:|---:|---:|
| QPLIB | 19 | 13 | 0.954 | 0.934 | 0.975 |
| Maros-Meszaros | 139 | 125 | 0.876 | 0.764 | 0.946 |
| Combined | 158 | 138 | 0.883 | 0.779 | 0.949 |

| Collection | Java valid inputs | Native valid inputs | Java only | Native only |
|---|---:|---:|---:|---:|
| QPLIB | 13 | 13 | 0 | 0 |
| Maros-Meszaros | 127 | 128 | 0 | 1 |
| Combined | 140 | 141 | 0 | 1 |

Coverage counts refer to the final recorded matched profile and include input aliases;
means count identical matrices only once. For every pair the native setup was additionally measured against the internal `info.setup_time` of OSQP; the difference is Python adapter work, so the matched setup + solve column removes it, and untouched solve ratios stay as measured.

The means describe the jointly validated solution set, not solver coverage on the excluded
inputs. The separately validated advanced Java profile for QPLIB 10038 above is not a
matched native profile; exclusion from this comparison does not undo that result.
A limit status leaves solvability unresolved. Even matching numerical infeasibility
verdicts do not constitute an independent mathematical proof. Excluded rows are screening
or diagnostic timings, and have no performance ratio.

#### QPLIB

| Instance | n | Java solve | Native solve | Java / native | Profile | Status / validation |
|---|---:|---:|---:|---:|---|---|
| `qplib_10034` | 40400 | 225.004839 s | 224.885370 s | — | A / 1e-09 | TIME_LIMIT_REACHED / RUN_TIME_LIMIT_REACHED; termination status |
| `qplib_10038` | 160800 | — | 224.064795 s | — | A / 1e-09 | RSS_MONITOR_UNAVAILABLE / RUN_TIME_LIMIT_REACHED; worker/resource failure |
| `qplib_8495` | 27543 | 0.467448 s | 0.476646 s | 0.981 | A | SOLVED / SOLVED |
| `qplib_8500` | 250997 | 225.016759 s | 221.303718 s | — | A | TIME_LIMIT_REACHED / RUN_TIME_LIMIT_REACHED; termination status |
| `qplib_8515` | 16002 | 225.001284 s | 224.995988 s | — | A | TIME_LIMIT_REACHED / RUN_TIME_LIMIT_REACHED; termination status |
| `qplib_8547` | 1003001 | — | 201.146007 s | — | A | RSS_MONITOR_UNAVAILABLE / RUN_TIME_LIMIT_REACHED; worker/resource failure |
| `qplib_8559` | 10000 | 32.475134 s | 35.983137 s | 0.903 | A | SOLVED / SOLVED |
| `qplib_8567` | 10000 | 53.632004 s | 59.984141 s | 0.894 | A | SOLVED / SOLVED |
| `qplib_8602` | 34552 | 31.716734 s | 32.483907 s | 0.976 | A | SOLVED / SOLVED |
| `qplib_8616` | 13870 | 0.051052 s | 0.050715 s | 1.007 | A | SOLVED / SOLVED |
| `qplib_8785` | 10399 | 2.008472 s | 1.981384 s | 1.014 | A | SOLVED / SOLVED |
| `qplib_8790` | 39204 | 0.429001 s | 0.441427 s | 0.972 | A | SOLVED / SOLVED |
| `qplib_8792` | 15129 | 0.377277 s | 0.373852 s | 1.009 | A | SOLVED / SOLVED |
| `qplib_8845` | 1546 | 7.551524 s | 8.292000 s | 0.911 | A | SOLVED / SOLVED |
| `qplib_8906` | 5223 | 0.966887 s | 1.066671 s | 0.906 | A | SOLVED / SOLVED |
| `qplib_8938` | 4001 | 56.896101 s | 64.030477 s | 0.889 | F10 | SOLVED / SOLVED |
| `qplib_8991` | 14400 | 0.188119 s | 0.185065 s | 1.016 | A | SOLVED / SOLVED |
| `qplib_9002` | 2890 | 9.517382 s | 10.183784 s | 0.935 | A | SOLVED / SOLVED |
| `qplib_9008` | 1009306 | — | — | — | A | PROCESS_ERROR / PROCESS_ERROR; worker/resource failure |

#### Maros-Meszaros

| Instance | n | Java solve | Native solve | Java / native | Profile | Status / validation |
|---|---:|---:|---:|---:|---|---|
| `aug2d` | 20200 | 0.027040 s | 0.030699 s | 0.881 | A | SOLVED / SOLVED |
| `aug2dc` | 20200 | 0.027107 s | 0.030475 s | 0.889 | A | SOLVED / SOLVED |
| `aug2dcqp` | 20200 | 0.932717 s | 1.023989 s | 0.911 | A | SOLVED / SOLVED |
| `aug2dqp` | 20200 | 0.838724 s | 0.924962 s | 0.907 | A | SOLVED / SOLVED |
| `aug3d` | 3873 | 0.004620 s | 0.004752 s | 0.972 | A | SOLVED / SOLVED |
| `aug3dc` | 3873 | 0.004589 s | 0.004832 s | 0.950 | A | SOLVED / SOLVED |
| `aug3dcqp` | 3873 | 0.010159 s | 0.010433 s | 0.974 | A | SOLVED / SOLVED |
| `aug3dqp` | 3873 | 0.010092 s | 0.010090 s | 1.000 | A | SOLVED / SOLVED |
| `boyd1` | 93261 | 225.007218 s | 224.954797 s | — | A | SOLVED_INACCURATE / SOLVED_INACCURATE; original-matrix tolerance failed |
| `boyd2` | 93263 | 225.006119 s | — | — | A | TIME_LIMIT_REACHED / RSS_MONITOR_UNAVAILABLE; worker/resource failure |
| `cont-050` | 2597 | 0.306124 s | 0.332666 s | 0.920 | A | SOLVED / SOLVED |
| `cont-100` | 10197 | 7.721537 s | 8.366774 s | 0.923 | A | SOLVED / SOLVED |
| `cont-101` | 10197 | 12.065546 s | 11.307979 s | 1.067 | A | SOLVED / SOLVED |
| `cont-200` | 40397 | 180.890199 s | 206.580491 s | 0.876 | A | SOLVED / SOLVED |
| `cont-201` | 40397 | 225.003876 s | 224.757087 s | — | A | TIME_LIMIT_REACHED / RUN_TIME_LIMIT_REACHED; termination status |
| `cont-300` | 90597 | 225.016819 s | 224.039014 s | — | A | TIME_LIMIT_REACHED / RUN_TIME_LIMIT_REACHED; termination status |
| `cvxqp1_l` | 10000 | 31.720245 s | 35.341094 s | 0.898 | A | SOLVED / SOLVED |
| `cvxqp1_m` | 1000 | 0.054985 s | 0.058324 s | 0.943 | A | SOLVED / SOLVED |
| `cvxqp1_s` | 100 | 0.002030 s | 0.001099 s | 1.848 | A | SOLVED / SOLVED |
| `cvxqp2_l` | 10000 | 4.894583 s | 5.390829 s | 0.908 | A | SOLVED / SOLVED |
| `cvxqp2_m` | 1000 | 0.013716 s | 0.013675 s | 1.003 | A | SOLVED / SOLVED |
| `cvxqp2_s` | 100 | 0.000168 s | 0.000178 s | 0.941 | A | SOLVED / SOLVED |
| `cvxqp3_l` | 10000 | 49.942269 s | 56.588915 s | 0.883 | A | SOLVED / SOLVED |
| `cvxqp3_m` | 1000 | 0.286179 s | 0.308428 s | 0.928 | A | SOLVED / SOLVED |
| `cvxqp3_s` | 100 | 0.000314 s | 0.000321 s | 0.978 | A | SOLVED / SOLVED |
| `dpklo1` | 133 | 0.000305 s | 0.000344 s | 0.886 | A | SOLVED / SOLVED |
| `dtoc3` | 14999 | 0.104994 s | 0.103936 s | 1.010 | A | SOLVED / SOLVED |
| `dual1` | 85 | 0.000235 s | 0.000264 s | 0.889 | A | SOLVED / SOLVED |
| `dual2` | 96 | 0.000285 s | 0.000328 s | 0.870 | A | SOLVED / SOLVED |
| `dual3` | 111 | 0.000369 s | 0.000429 s | 0.861 | A | SOLVED / SOLVED |
| `dual4` | 75 | 0.000188 s | 0.000212 s | 0.887 | A | SOLVED / SOLVED |
| `dualc1` | 9 | 0.001566 s | 0.001641 s | 0.954 | A | SOLVED / SOLVED |
| `dualc2` | 7 | 0.000596 s | 0.000594 s | 1.004 | A | SOLVED / SOLVED |
| `dualc5` | 8 | 0.000230 s | 0.000230 s | 0.997 | A | SOLVED / SOLVED |
| `dualc8` | 8 | 0.000392 s | 0.000402 s | 0.976 | A | SOLVED / SOLVED |
| `exdata` | 3000 | 1.874431 s | 1.982915 s | 0.945 | A | SOLVED / SOLVED |
| `genhs28` | 10 | 0.000010 s | 0.000018 s | 0.537 | A | SOLVED / SOLVED |
| `gouldqp2` | 699 | 0.028545 s | 0.027315 s | 1.045 | A | SOLVED / SOLVED |
| `gouldqp3` | 699 | 0.000781 s | 0.000678 s | 1.152 | A | SOLVED / SOLVED |
| `hs118` | 15 | 0.000201 s | 0.000245 s | 0.820 | A | SOLVED / SOLVED |
| `hs21` | 2 | 0.000004 s | 0.000013 s | 0.324 | A | SOLVED / SOLVED |
| `hs268` | 5 | 0.000038 s | 0.000044 s | 0.863 | A | SOLVED / SOLVED |
| `hs35` | 3 | 0.000004 s | 0.000011 s | 0.347 | A | SOLVED / SOLVED |
| `hs35mod` | 3 | 0.000004 s | 0.000011 s | 0.309 | A | SOLVED / SOLVED |
| `hs51` | 5 | 0.000006 s | 0.000013 s | 0.462 | A | SOLVED / SOLVED |
| `hs52` | 5 | 0.000006 s | 0.000012 s | 0.507 | A | SOLVED / SOLVED |
| `hs53` | 5 | 0.000003 s | 0.000010 s | 0.344 | A | SOLVED / SOLVED |
| `hs76` | 4 | 0.000005 s | 0.000012 s | 0.406 | A | SOLVED / SOLVED |
| `hues-mod` | 10000 | 0.210316 s | 0.327804 s | 0.642 | A | SOLVED / SOLVED |
| `huestis` | 10000 | 1.547928 s | 2.503274 s | 0.618 | A | SOLVED / SOLVED |
| `ksip` | 20 | 0.182049 s | 0.184860 s | 0.985 | A | SOLVED / SOLVED |
| `laser` | 1002 | 0.003141 s | 0.003110 s | 1.010 | A | SOLVED / SOLVED |
| `liswet1` | 10002 | 78.784634 s | 99.622786 s | 0.791 | A | SOLVED / SOLVED |
| `liswet10` | 10002 | 9.262303 s | 12.125113 s | 0.764 | A | SOLVED / SOLVED |
| `liswet11` | 10002 | 11.349060 s | 14.250865 s | 0.796 | A | SOLVED / SOLVED |
| `liswet12` | 10002 | 225.000834 s | 224.997350 s | — | A | TIME_LIMIT_REACHED / RUN_TIME_LIMIT_REACHED; termination status |
| `liswet2` | 10002 | 10.976067 s | 14.010017 s | 0.783 | A | SOLVED / SOLVED |
| `liswet3` | 10002 | 10.823428 s | 13.852229 s | 0.781 | A | SOLVED / SOLVED |
| `liswet4` | 10002 | 11.094500 s | 13.795301 s | 0.804 | A | SOLVED / SOLVED |
| `liswet5` | 10002 | 3.319049 s | 4.083575 s | 0.813 | A | SOLVED / SOLVED |
| `liswet6` | 10002 | 9.513291 s | 11.663250 s | 0.816 | A | SOLVED / SOLVED |
| `liswet7` | 10002 | 11.398909 s | 14.306348 s | 0.797 | A | SOLVED / SOLVED |
| `liswet8` | 10002 | 11.281033 s | 14.274626 s | 0.790 | A | SOLVED / SOLVED |
| `liswet9` | 10002 | 11.890721 s | 15.013930 s | 0.792 | A | SOLVED / SOLVED |
| `lotschd` | 12 | 0.000043 s | 0.000049 s | 0.891 | A | SOLVED / SOLVED |
| `mosarqp1` | 2500 | 0.027501 s | 0.028015 s | 0.982 | A | SOLVED / SOLVED |
| `mosarqp2` | 900 | 0.007312 s | 0.007977 s | 0.917 | A | SOLVED / SOLVED |
| `powell20` | 10000 | 119.070628 s | 110.427444 s | 1.078 | A | SOLVED / SOLVED |
| `primal1` | 325 | 0.000940 s | 0.001009 s | 0.931 | A | SOLVED / SOLVED |
| `primal2` | 649 | 0.000854 s | 0.000953 s | 0.895 | A | SOLVED / SOLVED |
| `primal3` | 745 | 0.002321 s | 0.002601 s | 0.892 | A | SOLVED / SOLVED |
| `primal4` | 1489 | 0.003684 s | 0.003812 s | 0.966 | A | SOLVED / SOLVED |
| `primalc1` | 230 | 0.087051 s | 0.080769 s | 1.078 | F0.1 / 1e-09 | SOLVED / SOLVED |
| `primalc2` | 231 | 0.428961 s | 0.562174 s | 0.763 | F0.1 / 1e-09 | SOLVED / SOLVED |
| `primalc5` | 287 | 0.089164 s | 0.108354 s | 0.823 | F10 / 1e-09 | SOLVED / SOLVED |
| `primalc8` | 520 | 0.419664 s | 0.500743 s | 0.838 | F0.1 / 1e-09 | SOLVED / SOLVED |
| `q25fv47` | 1571 | 8.149327 s | 9.071031 s | 0.898 | A | SOLVED / SOLVED |
| `qadlittl` | 97 | 0.004424 s | 0.004207 s | 1.052 | A | SOLVED / SOLVED |
| `qafiro` | 32 | 0.000121 s | 0.000129 s | 0.941 | A | SOLVED / SOLVED |
| `qbandm` | 472 | 0.267218 s | 0.249850 s | 1.070 | A | SOLVED / SOLVED |
| `qbeaconf` | 262 | 0.031069 s | 0.031530 s | 0.985 | A | SOLVED / SOLVED |
| `qbore3d` | 315 | 9.562453 s | 9.215581 s | — | A | MAX_ITER_REACHED / MAXIMUM_ITERATIONS_REACHED; termination status |
| `qbrandy` | 249 | 0.327545 s | 0.334419 s | 0.979 | A | SOLVED / SOLVED |
| `qcapri` | 353 | 13.443705 s | 13.355763 s | — | A | MAX_ITER_REACHED / MAXIMUM_ITERATIONS_REACHED; termination status |
| `qe226` | 282 | 2.080378 s | 2.204754 s | 0.944 | F0.1 | SOLVED / SOLVED |
| `qetamacr` | 688 | 7.771960 s | 5.174581 s | 1.502 | A | SOLVED / SOLVED |
| `qfffff80` | 854 | 6.515417 s | 7.085614 s | 0.920 | A | SOLVED / SOLVED |
| `qforplan` | 421 | 1.040480 s | 1.119775 s | 0.929 | F10 | SOLVED / SOLVED |
| `qgfrdxpn` | 1092 | 11.158702 s | 8.538262 s | 1.307 | A | SOLVED / SOLVED |
| `qgrow15` | 645 | 24.329887 s | 25.305472 s | 0.961 | F0.1 | SOLVED / SOLVED |
| `qgrow22` | 946 | 32.366284 s | 33.698725 s | 0.960 | F0.1 | SOLVED / SOLVED |
| `qgrow7` | 301 | 12.187601 s | 2.974565 s | — | A | MAX_ITER_REACHED / SOLVED; termination status |
| `qisrael` | 142 | 0.691608 s | 0.754795 s | 0.916 | A | SOLVED / SOLVED |
| `qpcblend` | 83 | 0.004202 s | 0.004462 s | 0.942 | A | SOLVED / SOLVED |
| `qpcboei1` | 384 | 0.205034 s | 0.205099 s | 1.000 | A | SOLVED / SOLVED |
| `qpcboei2` | 143 | 5.613146 s | 5.521021 s | — | A | MAX_ITER_REACHED / MAXIMUM_ITERATIONS_REACHED; termination status |
| `qpcstair` | 467 | 0.079042 s | 0.077724 s | 1.017 | A | SOLVED / SOLVED |
| `qpilotno` | 2172 | 64.620045 s | 68.424535 s | 0.944 | A | SOLVED / SOLVED |
| `qptest` | 2 | 0.000005 s | 0.000012 s | 0.398 | A | SOLVED / SOLVED |
| `qrecipe` | 180 | 0.001061 s | 0.001010 s | 1.051 | A | SOLVED / SOLVED |
| `qsc205` | 203 | 0.000969 s | 0.000922 s | 1.051 | A | SOLVED / SOLVED |
| `qscagr25` | 500 | 0.239169 s | 0.209449 s | 1.142 | A | SOLVED / SOLVED |
| `qscagr7` | 140 | 3.713410 s | 3.284716 s | — | A | MAX_ITER_REACHED / MAXIMUM_ITERATIONS_REACHED; termination status |
| `qscfxm1` | 457 | 1.218184 s | 1.565845 s | 0.778 | A | SOLVED / SOLVED |
| `qscfxm2` | 914 | 2.694437 s | 2.683378 s | 1.004 | A | SOLVED / SOLVED |
| `qscfxm3` | 1371 | 3.045642 s | 3.034866 s | 1.004 | A | SOLVED / SOLVED |
| `qscorpio` | 358 | 8.924327 s | 8.172260 s | 1.092 | F10 | SOLVED / SOLVED |
| `qscrs8` | 1169 | 3.749373 s | 3.550324 s | 1.056 | A | SOLVED / SOLVED |
| `qscsd1` | 760 | 0.019816 s | 0.018993 s | 1.043 | A | SOLVED / SOLVED |
| `qscsd6` | 1350 | 0.363485 s | 0.356097 s | 1.021 | A | SOLVED / SOLVED |
| `qscsd8` | 2750 | 0.077234 s | 0.073587 s | 1.050 | A | SOLVED / SOLVED |
| `qsctap1` | 480 | 1.022640 s | 1.003112 s | 1.019 | F0.1 | SOLVED / SOLVED |
| `qsctap2` | 1880 | 0.530640 s | 0.535530 s | 0.991 | F0.1 | SOLVED / SOLVED |
| `qsctap3` | 2480 | 0.334878 s | 0.338187 s | 0.990 | F0.1 | SOLVED / SOLVED |
| `qseba` | 1028 | 5.099360 s | 4.285953 s | 1.190 | A | SOLVED / SOLVED |
| `qshare1b` | 225 | 1.229500 s | 4.479737 s | — | A | SOLVED / SOLVED; objective disagreement |
| `qshare2b` | 79 | 0.458195 s | 0.461204 s | 0.993 | F0.1 | SOLVED / SOLVED |
| `qshell` | 1775 | 0.477738 s | 0.500762 s | 0.954 | A | SOLVED / SOLVED |
| `qship04l` | 2118 | 0.165723 s | 0.189591 s | 0.874 | A | SOLVED / SOLVED |
| `qship04s` | 1458 | 0.061056 s | 0.061571 s | 0.992 | A | SOLVED / SOLVED |
| `qship08l` | 4283 | 1.667317 s | 1.819666 s | 0.916 | A | SOLVED / SOLVED |
| `qship08s` | 2387 | 0.290702 s | 0.311463 s | 0.933 | A | SOLVED / SOLVED |
| `qship12l` | 5427 | 0.872539 s | 0.966791 s | 0.903 | A | SOLVED / SOLVED |
| `qship12s` | 2763 | 0.174376 s | 0.182660 s | 0.955 | A | SOLVED / SOLVED |
| `qsierra` | 2036 | 0.446106 s | 0.425081 s | 1.049 | A | SOLVED / SOLVED |
| `qstair` | 467 | 1.594997 s | 1.637672 s | 0.974 | A | SOLVED / SOLVED |
| `qstandat` | 1075 | 0.170879 s | 0.158469 s | 1.078 | A | SOLVED / SOLVED |
| `s268` | 5 | 0.000038 s | 0.000045 s | 0.844 | A | SOLVED / SOLVED; alias of hs268 |
| `selfqmin` | 21920 | 225.007431 s | 224.434100 s | — | A | TIME_LIMIT_REACHED / RUN_TIME_LIMIT_REACHED; termination status |
| `stadat1` | 2001 | 1.420074 s | 1.411056 s | 1.006 | F10 | SOLVED / SOLVED |
| `stadat2` | 2001 | 28.835707 s | 28.928301 s | 0.997 | F10 | SOLVED / SOLVED |
| `stadat3` | 4001 | 59.264812 s | 59.928050 s | 0.989 | F10 | SOLVED / SOLVED |
| `stcqp1` | 4097 | 0.012850 s | 0.012909 s | 0.995 | A | SOLVED / SOLVED |
| `stcqp2` | 4097 | 0.024794 s | 0.024707 s | 1.004 | A | SOLVED / SOLVED |
| `tame` | 2 | 0.000003 s | 0.000010 s | 0.270 | A | SOLVED / SOLVED |
| `ubh1` | 18009 | 8.249890 s | 8.703291 s | 0.948 | A | SOLVED / SOLVED |
| `values` | 202 | 0.005908 s | 0.006225 s | 0.949 | A | SOLVED / SOLVED |
| `yao` | 2002 | 50.975447 s | 60.532665 s | — | A | MAX_ITER_REACHED / MAXIMUM_ITERATIONS_REACHED; termination status |
| `zecevic2` | 2 | 0.000007 s | 0.000015 s | 0.447 | A | SOLVED / SOLVED |

Raw runs, matrix checksums, profiles, quality checks and worker failure reasons are
retained externally in `josqp-tools/results/full-release-bench.json` and
`full-release-summary.json`. The companion files are benchmark inputs, not bundled runtime
dependencies. No Python code is included in this Java repository.

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
