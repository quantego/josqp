#!/usr/bin/env python3
"""Summarize retained stability runs; numerical audits are not optimality certificates."""
import collections
import json
import math
import statistics
from pathlib import Path

RESULTS = Path(__file__).resolve().parent / 'results'
SIBLING = RESULTS.parents[2] / 'josqp_benchmarks'


def read(name):
    return json.loads((RESULTS / ('stability-' + name + '.json')).read_text())


def outcomes(data):
    return dict(collections.Counter(r.get('failure') or r['rows'][0]['status']
                                    for r in data['runs']))


def references():
    values = {}
    for name in ('MarosMeszaros', 'QPLIB'):
        for line in (SIBLING / 'Results' / ('report_' + name + '.txt')).read_text().splitlines():
            columns = line.split('|')
            if len(columns) < 3:
                continue
            fields = columns[2].split()
            if not fields:
                continue
            try:
                values[columns[1].strip() + '.qps'] = float(fields[0])
            except ValueError:
                pass
    return values


def summarize():
    for mode, count in (('portable', 47), ('vector', 50), ('baseline', 24), ('pre-workspace', 40)):
        assert 'OK (%d tests)' % count in (RESULTS / ('stability-' + mode + '-tests.txt')).read_text(), mode
    assert 'Tests run: 6,  Failures: 6' in (RESULTS / 'stability-regressions-before.txt').read_text()
    before, after, strict = (read(n) for n in
                             ('before-default', 'final-default', 'final-strict'))
    grouped = collections.defaultdict(list)
    for run in strict['runs']:
        assert not run.get('failure'), run
        assert len(run['rows']) == 1, run
        grouped[run['file']].append(run)
    refs = references()
    table = []
    for name, runs in grouped.items():
        assert len(runs) == 3 and {r['fork'] for r in runs} == {0, 1, 2}, name
        assert len({r['data_sha256'] for r in runs}) == 1, name
        rows = [r['rows'][0] for r in runs]
        fields = ('objective', 'iterations', 'rho_updates', 'x_hash', 'y_hash', 'status_polish')
        assert len({tuple(r[k] for k in fields) for r in rows}) == 1, name
        for row in rows:
            assert row['status'] == 'SOLVED' and row['finite'] == 'true', name
            for key in ('objective', 'independent_objective', 'primal_violation',
                        'stationarity', 'complementarity', 'dual_cone_error', 'gap_expression'):
                assert math.isfinite(float(row[key])), (name, key)
            assert math.isclose(float(row['objective']), float(row['independent_objective']),
                                rel_tol=1e-9, abs_tol=1e-9), name
        row = rows[0]
        value = float(row['objective'])
        table.append(dict(file=name, n=int(row['n']), m=int(row['m']),
                          objective=value, historical_rounded_gurobi=refs[name],
                          relative_difference_to_rounded_reference=abs(value-refs[name])/max(abs(refs[name]), 1e-30),
                          iterations=int(row['iterations']), rho_updates=int(row['rho_updates']),
                          status_polish=int(row['status_polish']),
                          setup_ms_median=statistics.median(float(r['setup_ms']) for r in rows),
                          solve_ms_median=statistics.median(float(r['solve_ms']) for r in rows),
                          solve_ms_range=[min(float(r['solve_ms']) for r in rows), max(float(r['solve_ms']) for r in rows)],
                          **{k: float(row[k]) for k in ('primal_violation', 'stationarity', 'complementarity',
                                                       'dual_cone_error', 'gap_expression')}))
    assert len(table) == 8
    assert {r['file'] for r in before['runs']} == {r['file'] for r in after['runs']}
    by_before = {r['file']: r for r in before['runs']}
    for run in after['runs']:
        if 'data_sha256' in run:
            assert run['data_sha256'] == by_before[run['file']]['data_sha256'], run['file']
    ubh_before = next(r['rows'][0] for r in read('tune-strict')['runs'] if r['file'] == 'ubh1.qps')
    ubh_fixed = next(r['rows'][0] for r in read('tune-strict-fixed')['runs'] if r['file'] == 'ubh1.qps')
    summary = dict(
        solver_commit=strict['solver_commit'], baseline_commit=before['solver_commit'],
        default_before=outcomes(before), default_after=outcomes(after),
        default_settings=dict(eps_abs=0.001, eps_rel=0.001, polish=False, solve_limit_seconds=5),
        strict_settings=dict(eps_abs=1e-6, eps_rel=1e-6, infeasibility_eps=1e-8, polish=True,
                             max_iter=100000, solve_limit_seconds=15, external_limit_seconds=50),
        timing='Three independent cold JVMs per case; no warmup; medians exclude setup and parsing.',
        reference='Rounded historical Gurobi table only; no fresh native solver run or exact reference certificate.',
        strict_cases=table,
        polish_example=dict(file='ubh1.qps', before_objective=float(ubh_before['objective']),
                            corrected_objective=float(ubh_fixed['objective']), historical_rounded_gurobi=refs['ubh1.qps'],
                            conclusion='Incorrect polish candidate rejected; reference accuracy still unresolved.'),
        unresolved_normalized_runs=[dict(r['rows'][0], file=r['file'])
                                    for r in read('tune-normalized')['runs']],
        limitations=['SOLVED uses global relative tolerances and does not imply small absolute row errors.',
                     'Original-unit gap expression is diagnostic, not a rigorous bound when residuals are nonzero.',
                     'New performance points change tolerances and iteration counts; they are not speedup comparisons.',
                     'Identical 32-bit vector hashes do not prove bit-identical full vectors.'])
    checks = dict(strict_runs=24, selected_cases=8, all_solved=True, all_finite=True,
                  independent_objective_matches=True, identical_reported_objectives_iterations_and_hashes_across_forks=True,
                  unchanged_input_hashes_before_after=True, missing_inputs=['xdata.qps'],
                  regression_tests_before=dict(total=6, failures=6),
                  junit_current=dict(portable=47, vector=50),
                  junit_archived=dict(baseline_qplib=24, pre_workspace=40),
                  optimality_certificates_verified=False)
    (RESULTS / 'stability-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (RESULTS / 'stability-verification.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps(dict(outcomes=summary['default_after'], strict_cases=table), indent=2))


if __name__ == '__main__':
    summarize()
