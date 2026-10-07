# jOSQP

jOSQP is a fork of the quadratic programming solver [OSQP](http://osqp.org) (Operator Splitting
Quadratic Program), written entirely in Java. It has no runtime dependencies beyond Java 8, and its
settings, defaults and termination criteria follow OSQP 1.0.

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
&& &y \leq 5
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
double[] upper = new double[n];
Model.Builder builder = Model.getBuilder();
builder.addVariables(n, lower, upper);
builder.linearObjective(new int[] {0, 1}, new double[] {1.0, 2.0});
builder.quadraticObjective(new int[] {0}, new int[] {0}, new double[] {0.5});
builder.constraint(new int[] {0, 1}, new double[] {2.0, -1.0}, 4.0, Double.POSITIVE_INFINITY);
Model model = builder.build();
```

For 300,000 variables, 2,000 rows and 200 entries per row this path builds the model in about 41 ms
and allocates 129 MB, against 224 ms and 236 MB along the algebraic builder.

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

`new OSQP(data, settings)` and `Model#solve(settings)` take an `OSQP.Settings` object that mirrors
`OSQPSettings` of OSQP 1.0. The defaults are the ones of the C library — `eps_abs` and `eps_rel`
`1e-3`, `max_iter` 4000, `rho` 0.1, `adaptive_rho_interval` 50, `check_dualgap` enabled,
`time_limit` `1e10`, polishing off, warm starting on — and an invalid configuration is rejected in
the `OSQP` constructor, as `osqp_setup` does.

Two families of names matter in practice:

* `polishing` and `warm_starting` are the OSQP 1.0 names and are available as accessors:
  `Settings#setPolishing(boolean)` / `isPolishing()` and `setWarmStarting(boolean)` /
  `isWarmStarting()`. The fields keep their 0.6 names `polish` and `warm_start`, so code written
  against older jOSQP keeps compiling.
* Settings that only exist in the C runtime are accepted for parity but cannot change a solve:
  `device`, `allocate_solution`, `profiler_level` and the four conjugate-gradient parameters
  (`cg_max_iter`, `cg_tol_reduction`, `cg_tol_fraction`, `cg_precond`), because jOSQP implements only
  the direct QDLDL factorisation. `linsys_solver = LinSysSolver.OSQP_INDIRECT_SOLVER` and
  `rho_is_vec = false` are rejected with an `IllegalArgumentException` rather than silently ignored,
  since honouring them is impossible here.

`OSQP.Info` carries the OSQP 1.0 fields `primdual_int` (running integral of the duality gap) and
`rel_kkt_error` alongside the usual status, objective and residual values.

Three jOSQP additions have no counterpart in OSQP 1.0 and default to the upstream behaviour:
`infeasibility_streak` and `infeasibility_margin` (see below), plus the optional diagnostics
`check_linear_residual` with `linear_refine_iter`, and `adaptive_rho_kkt`.

### Infeasibility certificates

OSQP decides primal and dual infeasibility from the last dual or primal iterate step: that step is
accepted as a certificate as soon as its residual passes `eps_prim_inf` / `eps_dual_inf`. On a
feasible problem a single step can cross that tolerance by numerical accident, and on badly scaled
problems the residual can plateau just below it for hundreds of checks. At the default tolerances
this makes both jOSQP and the upstream C solver report feasible instances as infeasible:

| Instance | Upstream verdict | jOSQP with the guards |
|---|---|---|
| `hues-mod` | `PRIMAL_INFEASIBLE` after 300 iterations | `SOLVED`, 875 iterations, objective `3.4761704979875020e7` |
| `huestis` | `PRIMAL_INFEASIBLE` after 100 iterations | `SOLVED`, 6675 iterations, objective `3.4755179890390e11` |
| `primalc5` | `DUAL_INFEASIBLE` | `MAX_ITER_REACHED` (as native OSQP at `eps_prim_inf = 1e-5`) |

Two settings guard the verdicts, and both have the upstream value available for exact parity:

- `Settings#infeasibility_streak` (default `2`): the certificate has to hold over that many
  consecutive termination checks, which rejects steps that only cross the tolerance once.
- `Settings#infeasibility_margin` (default `0.1`): safety factor on the effective certificate
  tolerance. Genuine certificates drive their residual to machine precision, so tightening the test
  by an order of magnitude keeps them while rejecting the plateaus. Set both to `1` to reproduce
  upstream OSQP verdicts exactly.

On all 63 benchmark instances the guards remove exactly the three false infeasibility verdicts above
and change no other status; the genuine infeasibility problems of the test suite are still detected
at the same iteration.

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

### The current version against native OSQP

Nine QPLIB instances through the same exported matrices, median of three cold-start runs, on an
Apple M4 with JDK 26. Both solvers use `eps_abs = eps_rel = 1e-3`, a 20 second limit and
`adaptive_rho_interval = 100`. jOSQP runs the plain Java 8 artifact (no Vector API backend) with
polishing off, native OSQP 1.1.3 with its default polishing, `eps_prim_inf = eps_dual_inf = 1e-8`
against 1e-4 in jOSQP, and 100,000 against 200,000 iterations.

| Instance | n | m | jOSQP | native OSQP | jOSQP status | native status |
|---|---:|---:|---:|---:|---|---|
| `qplib_8845` | 1,546 | 2,323 | 343 ms | 346 ms | SOLVED | solved |
| `qplib_9002` | 2,890 | 4,539 | 5,898 ms | 5,984 ms | SOLVED | solved |
| `qplib_8938` | 4,001 | 16,000 | over the limit | 14,747 ms | TIME_LIMIT_REACHED | maximum iterations reached |
| `qplib_8906` | 5,223 | 6,061 | 169 ms | 173 ms | SOLVED | solved |
| `qplib_8559` | 10,000 | 15,000 | 4,016 ms | 4,912 ms | SOLVED | solved |
| `qplib_8616` | 13,870 | 24,274 | 61 ms | 54 ms | SOLVED | solved |
| `qplib_8991` | 14,400 | 14,400 | 68 ms | 73 ms | SOLVED | solved |
| `qplib_8792` | 15,129 | 15,129 | 68 ms | 75 ms | SOLVED | solved |
| `qplib_8602` | 34,552 | 87,535 | over the limit | 12,993 ms | TIME_LIMIT_REACHED | run time limit reached, solved inaccurate |

On seven of the nine instances a pure Java solver stays with the C library, in five of them slightly
ahead of it. The two exceptions are honest gaps: `qplib_8602` is the largest case of the set and
native OSQP needs 13 seconds for it, so the 20 second limit is simply tight, and `qplib_8938` makes
both solvers give up. Reproduce with

```
../.tools/native-env/bin/python benchmarks/native_bench.py \
    --files ../josqp_benchmarks/extracted/qplib_8845.qps ... --solvers osqp \
    --repeats 3 --eps 1e-3 --time-limit 20 --output benchmarks/results/readme-native-osqp.json
java -cp target/portable-classes -Djosqp.bench.timeLimit=20 \
    com.quantego.josqp.SolverBenchmark ../josqp_benchmarks/extracted/<instance>.qps 0 3 200000
../.tools/native-env/bin/python benchmarks/summarize_readme.py
```

### What the current version changed

- The index-based bulk builder shortens model construction: 300,000 variables, 2,000 rows and 200
  entries per row build in about 41 ms and 129 MB instead of 224 ms and 236 MB.
- The optional Vector API backend accelerates the individual norm kernels by up to 7.5x, but its
  effect on a complete solve stays below one percent. The published artifact is therefore the plain
  Java 8 build.
- Reusing the symbolic factorisation and its work arrays across rho updates cut the solve time of
  `qplib_8602` from 31.6 to 27.2 seconds in an earlier measurement and removed almost all temporary
  allocations.

For the historical cross-solver comparison with Gurobi, OSQP and jOSQP across the whole
Maros-Meszaros and QPLIB sets see the [companion
repository](https://github.com/FaridAlavi/josqp_benchmarks); those numbers were measured on an
earlier version of jOSQP and are not comparable with the table above.

## Installation

The project is available at the central repository. Simply add the following dependency to your pom file
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
mvn test                            # 96 tests, JUnit 5
mvn -Pvector test                   # 99 tests, needs JDK 26 and jdk.incubator.vector
mvn -Prelease package               # main, sources and javadoc jar in target/
java -jar target/josqp-1.0.0.jar src/test/resources/qafiro.qps
```

`benchmarks/run.py` builds and tests the same sources without Maven and fetches the JUnit Platform
console launcher into the parent workspace on first use.

## Requirements

* Java 8 or higher to use the library; the published jar is compiled with `--release 8`
* JDK 26 for the optional Vector API profile
* A JDK and Maven 3.9 or newer to build and test

## Status

Version 1.0.0 fixes [issue #2](https://github.com/quantego/josqp/issues/2) (quadratic coefficients)
and [issue #4](https://github.com/quantego/josqp/issues/4) (updating `A` and `P`), and brings the
settings surface in line with OSQP 1.0 ([issue
#1](https://github.com/quantego/josqp/issues/1)) except the consolidated
`osqp_update_vec`/`osqp_update_data_mat` calls. [CHANGELOG.md](CHANGELOG.md) lists everything,
including the breaking changes.

## License

MIT License, see [LICENSE](LICENSE).

## Authors
* Nils Loehndorf
* Farid Alavi
