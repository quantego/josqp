# jOSQP

jOSQP is a fork of the quadratic programming solver, [OSQP](http://osqp.org) (Operator Splitting Quadratic Program), but written entirely in Java. There are no dependencies other than Java 8 or higher.

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
var builder = Model.getBuilder();
var x = builder.addVariable().lb(1).ub(3);
var y = builder.addVariable().ub(5);
builder.setObjective().add(5,x).add(0.5,x,x).add(3,y).minimize();
builder.addConstraint().add(2,x).add(-1,y).geq(4);
var model = builder.build();
model.solve()
```

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

### MPS/QPS Reader

Also, jOSQP comes with a parser for MPS and QPS files to import and solve existing linear and quadratic programming problems. Use the static function `Parser.readQmps` to read a file and pass the resulting model data to the solver:
```
var p = Parser.readQmps(<qps_file_name>);
var opt = new OSQP(p.getData(),new OSQP.Settings());
opt.solve();
```
The solver can also be called from the terminal
```
java -jar josqp.jar <qps_file_name>
```
After solving the problem through the terminal, jOSQP creates a log file,`josqp.log`, that contains a summary of the solution process. The terminal command has been used for benchmarking.

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

## Benchmarks

The following table summarizes the results of benchmarks with the Maros-Meszaros problem set as well as the QPlib with the original OSQP library written in C as well as the commercial GUROBI solver.

|                         | Gurobi | OSQP  | JOSQP |
| :---                    | :---:  | :---: | :---: |
| Maros-Meszaros  | 1.00   | 1.14  | 1.19  |
| QPLIB  | 1.02   | 1.20  | 1.00  |

The table reports the [shifted geometric mean](https://plato.asu.edu/ftp/shgeom.html) for each solver across all instances that could be solved by OSQP. Further details can be found [in a companion repository](https://github.com/FaridAlavi/josqp_benchmarks).

For the considered instances, jOSQP performs at par with the original OSQP code. This result is quite astonishing, as it demonstrates that a pure Java implementation of a solver for mathematical optimization is not necessarily slower than the same implementation written in C or C++, as long as such code does not make excessive use of objects and reuses memory as much as possible.

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

Building and testing locally requires a JDK and Maven; the tests run with `mvn test`, the optional
Vector API backend with `mvn -Pvector test` on JDK 26. See [CHANGELOG.md](CHANGELOG.md) for the
changes since 0.6.5.

## Requirements

* Java JDK 8

## Authors
* Nils Loehndorf
* Farid Alavi
