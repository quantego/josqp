# JOSQP

JOSQP is a fork of the solver [OSQP](osqp.org) (Operator Splitting Quadratic Program) that is written entirely in Java. There are no special dependencies for JOSQP and only Java Run Time in the target computer is enough to use this solver.

JOSQP solves the quadratic programs of the following form:
$$\min_x \quad \frac{1}{2}x^\mathrm{T}Px + q^\mathrm{T}x + c$$
$$\mathrm{s.t. \quad \quad \quad l \leq Ax \leq u}$$
where $x\in\mathbb{R}^n$ represents the decision variables and $P$ is a positive semi-definite matrix ($P\in S^n_+$).

The internal algorithm of the solver is based on the ADMM (Alternating Direction Method of Multipliers) and it is described [HERE](https://arxiv.org/abs/1711.08013).

## Usage

The JAR file `josqp.jar` and a Java Runtime Environment are the only requirements to solve QP problems using JOSQP. This solver is equipped with a QPS file reader. To solve a quadratic problem stored in a QPS file, use the following command in the terminal:
```
java -jar josqp.jar <qps_file_name>
```
After solving the QPS problem, JOSQP creates a log file with the default name `josqp.log`.

## Bug Reports

Please report any issue via the [GitHub issue tracker](https://github.com/quantego/josqp/issues). Any issue or suggestion is welcome.

## Numerical Benchmarks

A numerical benchmarking system is developed in [HERE](https://github.com/FaridAlavi/josqp_benchmarks).