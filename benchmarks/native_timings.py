#!/usr/bin/env python3
"""Median setup/solve timings for native OSQP and PIQP on JOSQP's exported matrices.

Each repetition builds a new solver on the same exported problem, so the numbers
measure setup plus solve without reuse. Reuse the adapters from native_compare.py.
"""
import argparse
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from native_compare import read_export, solve

FIELDS = ('status', 'finite', 'setup_ms', 'solve_ms', 'objective',
          'primal_violation', 'stationarity', 'gap_expression')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--eps', type=float, default=1e-6)
    parser.add_argument('--solvers', nargs='+', default=['osqp', 'piqp'])
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    result = {'arguments': vars(args), 'cases': []}
    for file in args.files:
        data = read_export(file)
        case = {'file': file, 'n': data['n'], 'm': data['m'],
                'export_sha256': data['export_sha256'], 'solvers': {}}
        for name in args.solvers:
            entry = {'runs': []}
            for _ in range(args.repeats):
                try:
                    row = solve(data, name, args.eps)
                    entry['runs'].append({k: row.get(k) for k in FIELDS})
                except Exception as exc:  # keep failures visible instead of aborting
                    entry['runs'].append({'failure': type(exc).__name__, 'detail': str(exc)})
            measured = [r for r in entry['runs'] if 'solve_ms' in r]
            if measured:
                entry['median_setup_ms'] = statistics.median(r['setup_ms'] for r in measured)
                entry['median_solve_ms'] = statistics.median(r['solve_ms'] for r in measured)
                entry['median_total_ms'] = statistics.median(
                    r['setup_ms'] + r['solve_ms'] for r in measured)
                entry['min_solve_ms'] = min(r['solve_ms'] for r in measured)
                entry['max_solve_ms'] = max(r['solve_ms'] for r in measured)
                entry['statuses'] = sorted({r['status'] for r in measured})
            case['solvers'][name] = entry
        result['cases'].append(case)

    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    compact = {c['file']: {s: {k: v for k, v in e.items() if k != 'runs'}
                           for s, e in c['solvers'].items()} for c in result['cases']}
    print(json.dumps(compact, indent=2))


if __name__ == '__main__':
    main()
