#!/usr/bin/env python3
"""Check equivalent complete runs and summarize paired factor-workspace measurements."""
import hashlib
import json
import math
import statistics
import struct
from pathlib import Path

root = Path(__file__).resolve().parent / 'results'
data = json.loads((root/'workspaces.json').read_text())
forks = data['arguments']['forks']; samples = data['arguments']['samples']
summary = {'baseline_commit':data['baseline_commit'],'solver_commit':data['solver_commit'],
           'cases':[],'failures':[r for r in data['runs'] if r.get('failure')]}
for file in sorted(data['data_sha256']):
    runs = [r for r in data['runs'] if r['file']==file]
    rows = [row for r in runs for row in r['rows']]
    complete = (len(runs)==2*forks and len({(r['mode'],r['fork']) for r in runs})==2*forks
                and all(len(r['rows'])==samples and not r.get('failure') for r in runs))
    signatures = {(r['status'],r['iterations'],r['objective'],r['primal_residual'],r['dual_residual'],
                   r['x_hash'],r['y_hash'],r['rho_updates']) for r in rows}
    hashes = {r.get('solution_sha256') for r in runs}
    bits_equal = complete and None not in hashes and len(hashes)==1
    solution_finite = False
    if bits_equal:
        r = runs[0]
        dump = root.parents[1]/'target'/('workspaces-'+Path(file).stem+'-'+r['mode']+'-'+str(r['fork'])+'.bin')
        blob=dump.read_bytes(); n,m=struct.unpack_from('>ii',blob)
        assert hashlib.sha256(blob).hexdigest()==r['solution_sha256']
        solution_finite = len(blob)==8+8*(n+m) and all(math.isfinite(v) for (v,) in struct.iter_unpack('>d',blob[8:]))
    modes={}
    for mode in ['pre-workspace','vector']:
        chosen=[r for r in runs if r['mode']==mode and len(r['rows'])==samples and not r.get('failure')]
        if not chosen:continue
        modes[mode]={key:statistics.median(statistics.median(float(row[key]) for row in r['rows']) for r in chosen)
                     for key in ['setup_ms','solve_ms','allocated_bytes','rho_updates']}
        modes[mode]['total_ms']=statistics.median(statistics.median(float(row['setup_ms'])+float(row['solve_ms']) for row in r['rows']) for r in chosen)
        modes[mode]['fork_solve_ms']=[statistics.median(float(row['solve_ms']) for row in r['rows']) for r in chosen]
        allocations = [int(row['allocated_bytes']) for r in chosen for row in r['rows']]
        modes[mode]['allocated_bytes_range']=[min(allocations),max(allocations)]
        modes[mode]['forks']=len(chosen)
    case={'file':file,'complete':complete,'metrics_exactly_equal':len(signatures)==1,'solution_bits_equal':bits_equal,
          'solution_finite':solution_finite,'status':sorted({r['status'] for r in rows}),'modes':modes}
    if complete:
        before=modes['pre-workspace'];after=modes['vector']
        case.update(solve_speedup=before['solve_ms']/after['solve_ms'],total_speedup=before['total_ms']/after['total_ms'],
                    solve_time_reduction_percent=100*(1-after['solve_ms']/before['solve_ms']))
    summary['cases'].append(case)
eligible=[c for c in summary['cases'] if c['complete'] and c['metrics_exactly_equal'] and c['solution_bits_equal']
          and c['solution_finite'] and c['status']==['SOLVED'] and c['file'].startswith('qplib_')]
summary['eligible_qplib_count']=len(eligible)
for metric in ['solve_speedup','total_speedup']:
    summary['qplib_geometric_mean_'+metric]=math.exp(statistics.mean(math.log(c[metric]) for c in eligible)) if eligible else None
(root/'workspaces-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
for c in summary['cases']:
    print(c['file'],'complete=',c['complete'],'bits=',c['solution_bits_equal'],
          {m:(round(v['solve_ms'],3),v['allocated_bytes']) for m,v in c['modes'].items()})
print({k:v for k,v in summary.items() if k.startswith('qplib_') or k=='eligible_qplib_count'})
