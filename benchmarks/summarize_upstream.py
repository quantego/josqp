#!/usr/bin/env python3
"""Check retained upstream experiments; aggregate per-JVM medians equally."""
import json
import math
from pathlib import Path
import re
import statistics
import struct
from collections import Counter

RESULTS = Path(__file__).resolve().parent / 'results'


def read(name):
    return json.loads((RESULTS / (name + '.json')).read_text())


def factor_summary(name, solver=False):
    data = read(name)
    summaries = []
    for file in data['arguments']['files']:
        runs = [r for r in data['runs'] if r['file'] == file]
        assert len(runs) == 4 * data['arguments']['forks']
        assert all(not r.get('failure') and r['exit_code'] == 0 and r['rows'] for r in runs)
        assert len({r['data_sha256'] for r in runs}) == 1
        hashes = ({r['solution_sha256'] for r in runs} if solver else
                  {v['factor_sha256'] for r in runs for v in r['rows']})
        assert len(hashes) == 1, (file, 'different binary factors or solutions')
        if solver:
            for r in runs:
                dump=RESULTS.parents[1]/'target'/('upstream-factor-'+Path(file).stem+'-'+r['mode']+'-'+str(r['fork'])+'.bin')
                binary=dump.read_bytes()
                n,m=struct.unpack_from('>ii',binary)
                assert (n,m)==(int(r['rows'][0]['n']),int(r['rows'][0]['m']))
                assert len(binary)==8+8*(n+m)
                assert all(math.isfinite(v[0]) for v in struct.iter_unpack('>d',binary[8:]))
            keys = ['status', 'iterations', 'objective', 'primal_residual',
                    'dual_residual', 'x_hash', 'y_hash', 'rho_updates']
            assert len({tuple(v[k] for k in keys) for r in runs for v in r['rows']}) == 1
            assert all(v['status'] == 'SOLVED' and
                       all(math.isfinite(float(v[k])) for k in keys[2:5])
                       for r in runs for v in r['rows'])
        timings = {}
        for mode in ['baseline', 'schedule', 'blocks', 'blocks-scalar']:
            selected = [r for r in runs if r['mode'] == mode]
            fields = ['solve_ms', 'setup_ms', 'allocated_bytes'] if solver else ['factor_ms']
            timings[mode] = {field: statistics.median([
                statistics.median(float(v[field]) for v in r['rows']) for r in selected
            ]) for field in fields}
            if solver:
                timings[mode]['setup_plus_solve_ms'] = statistics.median([
                    statistics.median(float(v['setup_ms']) + float(v['solve_ms'])
                                      for v in r['rows']) for r in selected])
        field = 'solve_ms' if solver else 'factor_ms'
        base = timings['baseline'][field]
        summaries.append(dict(file=file,identical_sha256=next(iter(hashes)),modes=timings,
                              time_saved_percent={m: 100*(1-v[field]/base)
                                                  for m,v in timings.items()}))
    return dict(cases=summaries,jvms=len(data['runs']),
                samples=sum(len(r['rows']) for r in data['runs']),
                geometric_time_saved_percent={m:100*(1-math.prod(
                    s['modes'][m]['solve_ms' if solver else 'factor_ms'] /
                    s['modes']['baseline']['solve_ms' if solver else 'factor_ms']
                    for s in summaries)**(1/len(summaries)))
                    for m in ['schedule','blocks','blocks-scalar']})


if __name__ == '__main__':
    summary = dict(factor=factor_summary('upstream-factor-micro'),
                   solver=factor_summary('upstream-factor-solvers',True))
    for name in ['upstream-gap-repeat','upstream-gap-disabled','upstream-linear-audit-final','upstream-kkt-rho',
                 'upstream-kkt-rho-strict','upstream-kkt-rho-long','upstream-kkt-rho-absolute',
                 'upstream-rho-control-absolute']:
        data=read(name)
        assert all(not r.get('failure') and r['exit_code']==0 and len(r['rows'])==1
                   for r in data['runs'])
        summary[name]={}
        for file in data['arguments']['files']:
            rows=[r['rows'][0] for r in data['runs'] if r['file']==file]
            assert len(rows)==data['arguments']['forks']
            assert len({tuple(v[k] for k in ['status','objective','iterations','rho_updates','x_hash','y_hash'])
                        for v in rows})==1, (name,file,'nonrepeatable outcomes')
            summary[name][file]=dict(statuses=[v['status'] for v in rows],
                median_solve_ms=statistics.median(float(v['solve_ms']) for v in rows),
                objective=[float(v['objective']) for v in rows],
                primal_violation=[float(v['primal_violation']) for v in rows],
                stationarity=[float(v['stationarity']) for v in rows],
                iterations=[int(v['iterations']) for v in rows],
                rho_updates=[int(v['rho_updates']) for v in rows])
    tests={}
    for mode in ['portable','vector']:
        log=(RESULTS/('upstream-final-'+mode+'-tests.txt')).read_text()
        match=re.search(r'OK \((\d+) tests\)',log)
        assert match,mode
        tests[mode]=int(match.group(1))
    native=read('upstream-native-status')
    exports=read('upstream-native-export')
    assert [(r['n'],r['m'],r['export_sha256']) for r in native['runs']]==[
        (r['n'],r['m'],r['export_sha256']) for r in exports['runs']]
    full=read('upstream-full-default')
    summary['full_default_statuses']=dict(Counter(
        r.get('failure') or r['rows'][0]['status'] for r in full['runs']))
    summary['full_default_records']=len(full['runs'])
    summary['native_adapters_executed']=native['adapters_executed']
    summary['tests']=tests
    (RESULTS/'upstream-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (RESULTS/'upstream-verification.json').write_text(json.dumps(dict(
        passed=True,tests=tests,full_vector_hashes_identical=True,
        full_vector_exports_finite=summary['solver']['jvms'],
        original_data_exports_checked=len(exports['runs']),
        factor_hashes_identical=True,native_adapters_executed=native['adapters_executed'],
        scope='Three QPLIB cases; SHA covers L values/D/Dinv or exported first-sample x/y.'),indent=2)+'\n')
    print(json.dumps(summary,indent=2))
