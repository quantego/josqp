# jOSQP

jOSQP is a fork of the quadratic programming solver, [OSQP](http://osqp.org) (Operator Splitting Quadratic Program), but written entirely in Java. There are no dependencies other than a JRE that is compatible with Java 8 or higher.

jOSQP solves quadratic programming problems of the form: 
$$\min_x \quad \frac{1}{2}x^\mathrm{T}Px + q^\mathrm{T}x + c$$
$$\mathrm{s.t. \quad \quad \quad l \leq Ax \leq u}$$
where $x\in\mathbb{R}^n$ represents variables and $P$ is a positive semi-definite matrix ($P\in S^n_+$).

The internal algorithm of the solver is based on the ADMM (Alternating Direction Method of Multipliers) and it is described [in the original OSQP paper](https://arxiv.org/abs/1711.08013).

## Usage

Unlike OSQP, jOSQP comes with an algebraic model builder that enables expressing models in a natural way.

### Model Builder

To make model building as easy as possible, jOSQP offers to model variables, constraints, and linear expressions as Java objects. Unfortunately, there is no operator overloading in Java, so that use of mathematical operators is not ossible.

Nonetheless, decision variables as well as their bounds, constraint, and objective function can be created in a single line. For example, a simple quadratic program such as

$$
\min 5x + 0.5x^2 + 3y \\
s.t. 2x-y \geq 4 \\
1 \leq  x \leq 3 \\
y \leq 5 \\
$$

can be expressed as

```
var builder = Model.getBuilder();
var x = builder.addVariable().lb(1).ub(3);
var y = builder.addVariable().ub(5);
builder.setObjective().add(5,x).add(0.5,x,x).add(3,y);
builder.addConstraint().add(2,x).add(-1,y).leq(3);
var model = builder.build();
```


### QPS Reader

The solver is comes with a parser for MPS and QPS files to import and solve existing problems. Use the static function `Parser.readQmps` to read a file and pass the resulting model data to the solver.
```
var p = Parser.readQmps("/home/user/myqps.qps");
var opt = new OSQP(p.getData(),new OSQP.Settings());
opt.solve();
```
The solver can also be called from the terminal for benchmarking using
```
java -jar josqp.jar <qps_file_name>
```
After solving the problem, jOSQP creates a log file with the default name `josqp.log`.

## Benchmarks

Extensive benchmarks with the original OSQP library written in C as well as the commercial GUROBI solver can be found under this link [here](https://github.com/FaridAlavi/josqp_benchmarks). For a sample of instances of the QPLIB, jOSQP performs at par with the original OSQP code, which demonstrates that a pure Java implementation of a solver for mathematical optimization is not necessarily slower than a native
implementation written in C or C++.
## Installation

The project is available at the central repository. Simply add the following dependency to your pom file
```
<dependency>
  <groupId>com.quantego</groupId>
  <artifactId>josqp</artifactId>
  <version>0.6.0</version>
</dependency>
```

Or you can download the [latest build](https://github.com/quantego/clp-java/releases/latest) of the jar from the release page and import it in your Java project.


## Requirements

* Java JDK 8

## Authors
* Nils Loehndorf
* Farid Alavi
