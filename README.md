# jOSQP

jOSQP is a fork of the quadratic programming solver, [OSQP](http://osqp.org) (Operator Splitting Quadratic Program), but written entirely in Java. There are no dependencies other than a JRE that is compatible with Java 8 or higher.

jOSQP solves quadratic programming problems of the form: 
$$\min_x \quad \frac{1}{2}x^\mathrm{T}Px + q^\mathrm{T}x + c$$
$$\mathrm{s.t. \quad \quad \quad l \leq Ax \leq u}$$
where $x\in\mathbb{R}^n$ represents variables and $P$ is a positive semi-definite matrix ($P\in S^n_+$).

The internal algorithm of the solver is based on the ADMM (Alternating Direction Method of Multipliers) and it is described [in the original OSQP paper](https://arxiv.org/abs/1711.08013).

## Usage

While OSQP has been linked with a few Python packages to facilitate model creation, jOSQP provides its own algebraic model builder as well as a parser for MPS and QPS files.

### Model Builder

To make model building as easy as possible, jOSQP offers to model variables, constraints, and linear expressions as Java objects. Unfortunately, there is no operator overloading in Java, so that use of mathematical operators is not possible. Nonetheless, decision variables as well as their bounds, constraint, and objective function can be created in a single line. For example, a simple quadratic program such as

$$
\begin{align}
&\min\ &5x + 0.5x^2 + 3y \\
&s.t.\  &2x-y \geq 4 \\
&&&1 \leq  x \leq 3 \\
&&y \leq 5
\end{align}
$$

can be expressed by only a few lines of code.

```
var builder = Model.getBuilder();
var x = builder.addVariable().lb(1).ub(3);
var y = builder.addVariable().ub(5);
builder.setObjective().add(5,x).add(0.5,x,x).add(3,y);
builder.addConstraint().add(2,x).add(-1,y).leq(3);
var model = builder.build();
```


### MPS/QPS Reader

Additionally, jOSQP comes with a parser for MPS and QPS files to import and solve existing linear and quadratic programming problems. Use the static function `Parser.readQmps` to read a file and pass the resulting model data to the solver:
```
var p = Parser.readQmps(<qps_file_name>);
var opt = new OSQP(p.getData(),new OSQP.Settings());
opt.solve();
```
The solver can also be called from the terminal for benchmarking using
```
java -jar josqp.jar <qps_file_name>
```
After solving the problem through the terminal, jOSQP creates a log file,`josqp.log`, that contains a summary of the solution process.

## Benchmarks

The following table summarizes the results of benchmarks with the Maros-Meszaros problem set as well as the QPlib with the original OSQP library written in C as well as the commercial GUROBI solver.

|                         | Gurobi | OSQP  | JOSQP |
| :---                    | :---:  | :---: | :---: |
| Maros-Meszaros  | 1.00   | 1.14  | 1.19  |
| QPLIB  | 1.02   | 1.20  | 1.00  |

The table shows the [shifted geometric mean](https://plato.asu.edu/ftp/shgeom.html) for each solver across all instances. Further details can be found [in a companion repository](https://github.com/FaridAlavi/josqp_benchmarks).

For the considered instances, jOSQP performs at par with the original OSQP code. This result is quite astonishing, as it demonstrates that a pure Java implementation of a solver for mathematical optimization is not necessarily slower than the same implementation written in C or C++, as long as such code does not make excessive use of objects and reuses memory as much as possible.

## Installation

The project is available at the central repository. Simply add the following dependency to your pom file
```
<dependency>
  <groupId>com.quantego</groupId>
  <artifactId>josqp</artifactId>
  <version>0.6.5</version>
</dependency>
```

Or you can download the [latest build](https://github.com/loehndorf/josqp/releases/latest) of the jar from the release page and import it in your Java project.

## Requirements

* Java JDK 8

## Authors
* Nils Loehndorf
* Farid Alavi
