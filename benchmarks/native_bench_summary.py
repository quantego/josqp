#!/usr/bin/env python3
"""Compare native OSQP with JOSQP on the full selected case list.

Both sides use their own solver time (without parsing). JOSQP is expected to be
measured with the strict profile (eps 1e-6, polishing, gap check) and the native
adapter with the same tolerances. Only cases both solvers actually converge to
the same objective enter the geometric mean of the time ratio.
"""
import argparse
import json
import math
import statistics
from pathlib import Path

JOSQP_CONVERGED = {'SOLVED', 'SOLVED_INACCURATE'}
OSQP_CONVERGED = {'solved', 'solved inaccurate'}


def median(values):
    return statistics.median(values) if values else None


def load_native(path):
    data = json.loads(Path(path).read_text())
    cases = {}
    for case in data['cases']:
        entry = case['solvers'].get('osqp', {})
        runs = [r for r in entry.get('runs', []) if 'solve_ms' in r]
        if not runs:
            cases[case['name']] = {'statuses': entry.get('runs'), 'missing': True}
            continue
        objectives = [r['objective'] for r in runs if r.get('objective') is not None]
        cases[case['name']] = {
            'n': case['n'], 'm': case['m'],
            'solve_ms': median([r['solve_ms'] for r in runs]),
            'setup_ms': median([r['setup_ms'] for r in runs]),
            'statuses': sorted({r['status'] for r in runs}),
            'primal_violation': max(r.get('primal_violation') or 0.0 for r in runs),
            'objective': median(objectives) if objectives else None,
        }
    return cases


def load_josqp(path):
    data = json.loads(Path(path).read_text())
    cases = {}
    for run in data['runs']:
        name = run['file'][:-4] if run['file'].endswith('.qps') else run['file']
        rows = run.get('rows') or []
        if not rows:
            continue
        entry = cases.setdefault(name, {'rows': []})
        entry['rows'].extend(rows)
    for name, entry in cases.items():
        rows = entry['rows']
        entry['solve_ms'] = median([float(r['solve_ms']) for r in rows if r.get('solve_ms') is not None])
        entry['setup_ms'] = median([float(r['setup_ms']) for r in rows if r.get('setup_ms') is not None])
        entry['statuses'] = sorted({r['status'] for r in rows if r.get('status')})
        entry['objective'] = median([float(r['independent_objective']) for r in rows
                                     if r.get('independent_objective') is not None])
        entry['primal_violation'] = max([float(r['primal_violation']) for r in rows
                                         if r.get('primal_violation') is not None] or [0.0])
    return cases


def agrees(a, b, tolerance=1e-5):
    if a is None or b is None:
        return False
    return abs(a - b) <= tolerance * max(1.0, abs(a), abs(b))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--native', default='benchmarks/results/upstream-native-bench.json')
    parser.add_argument('--josqp', action='append',
                        default=None,
                        help='may be repeated; later files override earlier ones per case')
    parser.add_argument('--output', default='benchmarks/results/upstream-native-bench-summary.json')
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    native = load_native(root / args.native)
    files = args.josqp or ['benchmarks/results/upstream-native-josqp-strict-full.json',
                           'benchmarks/results/upstream-native-josqp-warm.json']
    josqp = {}
    for path in files:
        if not (root / path).exists():
            continue
        josqp.update(load_josqp(root / path))

    rows = []
    for name in sorted(set(native) | set(josqp)):
        n = native.get(name, {})
        j = josqp.get(name, {})
        row = {'name': name, 'n': n.get('n', j.get('n')), 'm': n.get('m', j.get('m')),
               'osqp_ms': n.get('solve_ms'), 'josqp_ms': j.get('solve_ms'),
               'osqp_setup_ms': n.get('setup_ms'), 'josqp_setup_ms': j.get('setup_ms'),
               'osqp_status': ','.join(n.get('statuses') or []),
               'josqp_status': ','.join(j.get('statuses') or []),
               'osqp_objective': n.get('objective'), 'josqp_objective': j.get('objective')}
        osqp_statuses = set(n.get('statuses') or [])
        josqp_statuses = set(j.get('statuses') or [])
        row['osqp_converged'] = bool(osqp_statuses) and osqp_statuses <= OSQP_CONVERGED
        row['josqp_converged'] = bool(josqp_statuses) and josqp_statuses <= JOSQP_CONVERGED
        row['objective_agrees'] = agrees(row['osqp_objective'], row['josqp_objective'])
        row['comparable'] = bool(row['osqp_converged'] and row['josqp_converged']
                                 and row['objective_agrees'] and row['osqp_ms'] and row['josqp_ms'])
        if row['comparable']:
            row['ratio'] = row['josqp_ms'] / row['osqp_ms']
            row['ratio_total'] = ((row['josqp_ms'] + (row['josqp_setup_ms'] or 0.0))
                                  / (row['osqp_ms'] + (row['osqp_setup_ms'] or 0.0)))
        rows.append(row)

    comparable = [r for r in rows if r['comparable']]
    ratios = [r['ratio'] for r in comparable]
    totals = [r['ratio_total'] for r in comparable]

    summary = {
        'cases_measured': len(rows),
        'comparable_cases': len(comparable),
        'josqp_converged_cases': sum(1 for r in rows if r['josqp_converged']),
        'osqp_converged_cases': sum(1 for r in rows if r['osqp_converged']),
        'both_converged_but_objective_differs': sum(
            1 for r in rows if r['josqp_converged'] and r['osqp_converged'] and not r['objective_agrees']),
        'only_osqp_converged': sorted(r['name'] for r in rows if r['osqp_converged'] and not r['josqp_converged']),
        'only_josqp_converged': sorted(r['name'] for r in rows if r['josqp_converged'] and not r['osqp_converged']),
        'geometric_mean_ratio_solve': math.exp(sum(math.log(r) for r in ratios) / len(ratios)) if ratios else None,
        'geometric_mean_ratio_solve_plus_setup': math.exp(sum(math.log(r) for r in totals) / len(totals)) if totals else None,
        'median_ratio_solve': statistics.median(ratios) if ratios else None,
        'min_ratio_solve': min(ratios) if ratios else None,
        'max_ratio_solve': max(ratios) if ratios else None,
        'josqp_slower_count': sum(1 for r in ratios if r > 1.0),
        'josqp_faster_count': sum(1 for r in ratios if r <= 1.0),
        'geometric_std_solve': math.exp(statistics.pstdev([math.log(r) for r in ratios])) if len(ratios) > 1 else None,
        'by_native_time_threshold_ms': {
            str(threshold): {
                'cases': len(subset),
                'geometric_mean_ratio_solve': math.exp(sum(math.log(r['ratio']) for r in subset) / len(subset)),
                'median_ratio_solve': statistics.median([r['ratio'] for r in subset]),
                'josqp_faster_count': sum(1 for r in subset if r['ratio'] <= 1.0),
            }
            for threshold in (0, 20, 100, 500, 1000, 3000)
            if (subset := [r for r in comparable if r['osqp_ms'] >= threshold])
        },
    }
    payload = {'summary': summary, 'cases': rows}
    (root / args.output).write_text(json.dumps(payload, indent=2) + '\n')

    print(json.dumps(summary, indent=2))
    print('\n%-24s %10s %10s %8s  %-22s %-22s' % ('case', 'osqp ms', 'josqp ms', 'ratio', 'osqp', 'josqp'))
    for r in sorted(comparable, key=lambda x: -x['ratio']):
        print('%-24s %10.3f %10.3f %8.2f  %-22s %-22s' % (
            r['name'], r['osqp_ms'], r['josqp_ms'], r['ratio'], r['osqp_status'], r['josqp_status']))
    print('\nnot comparable:')
    for r in rows:
        if not r['comparable']:
            print('  %-24s osqp=%-24s josqp=%-24s' % (r['name'], r['osqp_status'], r['josqp_status']))


if __name__ == '__main__':
    main()
