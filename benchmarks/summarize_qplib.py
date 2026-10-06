#!/usr/bin/env python3
"""Summarize full forks only, separating incomplete solves and exact solution checks."""
import json
import math
import statistics
import struct
from pathlib import Path

root = Path(__file__).resolve().parent / 'results'
reference = root.parents[2] / 'josqp_benchmarks/Results/report_QPLIB.txt'
flagged = set()
for line in reference.read_text().splitlines():
    fields = line.split('|')
    if line.startswith('*') and len(fields) > 1 and fields[1].strip().startswith('qplib_'):
        flagged.add(fields[1].strip() + '.qps')
data = json.loads((root / 'qplib-compare.json').read_text())
files = sorted({r['file'] for r in data['runs']})
summary = {'cases': [], 'failures': [r for r in data['runs'] if r.get('failure')]}
for file in files:
    runs = [r for r in data['runs'] if r['file'] == file]
    rows = [row for r in runs for row in r['rows']]
    complete = len(runs) == 9 and all(len(r['rows']) == 2 and not r.get('failure') for r in runs)
    signatures = {(r['status'], r['iterations'], r['objective'], r['primal_residual'],
                   r['dual_residual'], r['x_hash'], r['y_hash']) for r in rows}
    hashes = {r.get('solution_sha256') for r in runs}
    finite = all(math.isfinite(float(r[k])) for r in rows
                 for k in ['objective', 'primal_residual', 'dual_residual'])
    solution_finite = False
    if complete and None not in hashes and len(hashes) == 1:
        r = runs[0]
        dump = root.parents[1] / 'target' / ('qplib-' + Path(file).stem + '-' + r['mode'] + '-' + str(r['fork']) + '.bin')
        blob = dump.read_bytes()
        n, m = struct.unpack_from('>ii', blob)
        solution_finite = len(blob) == 8 + 8 * (n + m) and all(math.isfinite(v) for (v,) in struct.iter_unpack('>d', blob[8:]))
    modes = {}
    for mode in ['baseline-qplib', 'scalar', 'vector']:
        chosen = [r for r in runs if r['mode'] == mode and len(r['rows']) == 2 and not r.get('failure')]
        if not chosen: continue
        medians = {key: [statistics.median(float(row[key]) for row in r['rows']) for r in chosen]
                   for key in ['setup_ms', 'solve_ms', 'allocated_bytes']}
        modes[mode] = {key: statistics.median(values) for key, values in medians.items()}
        modes[mode]['setup_plus_solve_ms'] = statistics.median(
            statistics.median(float(row['setup_ms']) + float(row['solve_ms']) for row in r['rows']) for r in chosen)
        modes[mode]['fork_solve_ms'] = medians['solve_ms']
        modes[mode]['forks'] = len(chosen)
    case = {'file': file, 'complete': complete, 'metrics_exactly_equal': len(signatures) == 1,
            'historical_reference_flag': file in flagged,
            'solution_bits_equal': complete and None not in hashes and len(hashes) == 1,
            'finite': finite, 'solution_finite': solution_finite,
            'status': sorted({r['status'] for r in rows}), 'modes': modes}
    if rows:
        case.update({k: rows[0][k] for k in ['n', 'm', 'nnz', 'iterations', 'objective', 'primal_residual', 'dual_residual']})
    if complete and len(modes) == 3:
        base = modes['baseline-qplib']['solve_ms']; scalar = modes['scalar']['solve_ms']; vector = modes['vector']['solve_ms']
        case.update(scalar_speedup=base/scalar, vector_speedup=base/vector,
                    extra_vector_speedup=scalar/vector, vector_time_reduction_percent=100*(1-vector/base))
    summary['cases'].append(case)
eligible = [c for c in summary['cases'] if c['complete'] and c['solution_bits_equal']
            and c['metrics_exactly_equal'] and c['finite'] and c['solution_finite'] and c['status'] == ['SOLVED']]
summary['eligible_solved_count'] = len(eligible)
for key in ['scalar_speedup', 'vector_speedup', 'extra_vector_speedup']:
    summary['geometric_mean_'+key] = math.exp(statistics.mean(math.log(c[key]) for c in eligible)) if eligible else None
primary = [c for c in eligible if not c['historical_reference_flag']]
summary['historically_unflagged_solved_count'] = len(primary)
for key in ['scalar_speedup', 'vector_speedup', 'extra_vector_speedup']:
    summary['unflagged_geometric_mean_'+key] = math.exp(statistics.mean(math.log(c[key]) for c in primary)) if primary else None
for group, label in [(eligible, 'all_solved'), (primary, 'unflagged')]:
    for mode in ['scalar', 'vector']:
        summary[label + '_geometric_mean_setup_plus_solve_' + mode + '_speedup'] = (
            math.exp(statistics.mean(math.log(c['modes']['baseline-qplib']['setup_plus_solve_ms'] /
                                              c['modes'][mode]['setup_plus_solve_ms']) for c in group)) if group else None)
(root / 'qplib-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
for c in summary['cases']:
    print(c['file'], c['status'], 'exact=', c['solution_bits_equal'],
          {m: round(v['solve_ms'], 3) for m, v in c['modes'].items()},
          'speedup=', round(c.get('vector_speedup', 0), 3))
print('SOLVED with identical outputs:', len(eligible))
print({k:v for k,v in summary.items() if k.startswith('geometric')})
