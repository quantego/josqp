#!/usr/bin/env python3
"""Build the benchmark table of the README from the two fresh runs.

  ../.tools/native-env/bin/python benchmarks/summarize_readme.py

Reads ``benchmarks/results/readme-native-osqp.json`` (native OSQP on the exported
matrices) and ``benchmarks/results/readme-josqp.csv`` (the same instances through
SolverBenchmark) and prints markdown with the median solve time per instance.
Both files are produced with the same tolerance (eps 1e-3) and a 20 second limit.
"""
import csv
import json
import statistics
from pathlib import Path

RESULTS = Path(__file__).resolve().parent / 'results'
NATIVE = RESULTS / 'readme-native-osqp.json'
JOSQP = RESULTS / 'readme-josqp.csv'


def native_by_case():
    data = json.loads(NATIVE.read_text())
    return {case['name']: case for case in data['cases']}, data['arguments']


def josqp_by_case():
    samples = {}
    with JOSQP.open() as handle:
        for row in csv.DictReader(handle):
            name = Path(row['file']).stem
            samples.setdefault(name, []).append(row)
    return samples


def main():
    native, arguments = native_by_case()
    josqp = josqp_by_case()
    print('| Instance | n | m | jOSQP solve | native OSQP solve | jOSQP status | native status |')
    print('|---|---:|---:|---:|---:|---|---|')
    for name in sorted(native, key=lambda n: native[n]['n']):
        case = native[name]
        entry = case['solvers']['osqp']
        runs = josqp.get(name, [])
        solved = [float(r['solve_ms']) for r in runs if r['status'] == 'SOLVED']
        josqp_ms = f'{statistics.median(solved):,.0f} ms' if solved else 'n/a'
        statuses = sorted({r['status'] for r in runs}) or ['n/a']
        print(f"| `{name}` | {case['n']:,} | {case['m']:,} | {josqp_ms} | "
              f"{entry['median_solve_ms']:,.0f} ms | {', '.join(statuses)} | "
              f"{', '.join(sorted(set(entry['statuses'])))} |")
    print()
    print(f"native arguments: eps={arguments['eps']}, repeats={arguments['repeats']}, "
          f"time limit={arguments['time_limit']}s")


if __name__ == '__main__':
    main()
