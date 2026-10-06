#!/usr/bin/env python3
"""Summarize by fork medians and verify full solver output matches across modes."""
import collections
import csv
import json
import pathlib
import statistics

ROOT = pathlib.Path(__file__).resolve().parent

def solvers(paths):
    cases = collections.defaultdict(list)
    for path in paths:
        for row in csv.DictReader(path.open()):
            cases[row['file']].append(row)
    result = {}
    for file, rows in cases.items():
        modes = collections.defaultdict(list)
        for row in rows:
            modes[row['mode']].append(row)
        entry = {'n': int(rows[0]['n']), 'm': int(rows[0]['m']), 'nnz': int(rows[0]['nnz']), 'modes': {}}
        for mode, samples in modes.items():
            forks = collections.defaultdict(list)
            for row in samples:
                forks[row['fork']].append(row)
            stats = {}
            for metric in ['setup_ms', 'solve_ms', 'total_ms', 'parse_ms', 'allocated_bytes']:
                def value(row):
                    return float(row['setup_ms']) + float(row['solve_ms']) if metric == 'total_ms' else float(row[metric])
                medians = [statistics.median(value(x) for x in fork) for fork in forks.values()]
                stats[metric] = statistics.median(medians)
                stats[metric + '_forks'] = medians
            stats.update(samples=len(samples), forks=len(forks), status=sorted({x['status'] for x in samples}),
                         iterations=sorted({int(x['iterations']) for x in samples}))
            entry['modes'][mode] = stats
        for metric in ['status', 'iterations', 'objective', 'primal_residual', 'dual_residual', 'x_hash', 'y_hash']:
            values = {row[metric] for row in rows}
            if len(values) != 1:
                raise AssertionError((file, metric, values))
        entry['outputs_match'] = True
        result[file] = entry
    return result

def kernels(paths):
    groups = collections.defaultdict(list)
    for path in paths:
        for row in csv.DictReader(path.open()):
            groups[(row['kernel'], row['n'], row['entries_per_column'])].append(row)
    return [{'kernel': key[0], 'n': int(key[1]), 'entries_per_column': int(key[2]),
             'scalar_ns': statistics.median(float(x['scalar_ns']) for x in rows),
             'vector_ns': statistics.median(float(x['vector_ns']) for x in rows),
             'speedup': statistics.median(float(x['speedup']) for x in rows), 'samples': len(rows)}
            for key, rows in sorted(groups.items())]

if __name__ == '__main__':
    results = ROOT / 'results'
    summary = {'solvers': solvers([results / 'solvers-hotloops.csv', results / 'qap15-hotloops.csv']),
               'before_triangular_optimization': solvers([results / 'solvers.csv', results / 'qap15.csv']),
               'kernels': kernels(sorted(results.glob('kernels-final-*.csv')))}
    (results / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary['solvers'], indent=2))
