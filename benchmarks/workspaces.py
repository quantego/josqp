#!/usr/bin/env python3
"""Paired comparison of the preceding optimized solver and reusable LDL workspaces."""
import argparse
import csv
import hashlib
import io
import json
import random
import subprocess
import time
from run import ROOT, JAVA, PRE_WORKSPACE, build

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', nargs='+')
    parser.add_argument('--forks', type=int, default=3)
    parser.add_argument('--samples', type=int, default=2)
    parser.add_argument('--warmup', type=int, default=1)
    parser.add_argument('--timeout', type=int, default=180)
    parser.add_argument('--output', default='benchmarks/results/workspaces.json')
    args = parser.parse_args()
    results = ROOT / 'benchmarks/results'
    if args.files is None:
        previous = json.loads((results/'qplib-summary.json').read_text())
        names = [c['file'] for c in previous['cases'] if c['complete'] and not c['historical_reference_flag'] and c['status']==['SOLVED']]
        names.append('qap15.mps')
    else:
        names = args.files
    files = []
    for name in names:
        file = ROOT.parent / 'josqp_benchmarks/extracted' / name
        if not file.exists(): file = ROOT / 'src/test/resources' / name
        if not file.exists(): raise FileNotFoundError(file)
        files.append(file)
    modes = ['pre-workspace', 'vector']
    for mode in modes: build(mode)
    data = {'arguments': vars(args), 'baseline_commit': PRE_WORKSPACE,
            'solver_commit': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'both_backends': 'vector', 'heap': '2g', 'max_iter':100000,
            'java':subprocess.run([str(JAVA/'bin/java'),'-version'],capture_output=True,text=True).stderr,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                             for directory in ['src/main/java','src/vector/java','benchmarks/java']
                             for p in sorted((ROOT/directory).rglob('*.java'))},
            'data_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files},'runs':[]}
    rng = random.Random(20261006)
    output = ROOT / args.output
    for file in files:
        for fork in range(args.forks):
            order = list(modes); rng.shuffle(order)
            for mode in order:
                dump = ROOT/'target'/('workspaces-'+file.stem+'-'+mode+'-'+str(fork)+'.bin')
                command = [str(JAVA/'bin/java'),'-Xms256m','-Xmx2g','--add-modules','jdk.incubator.vector',
                           '-Djosqp.vector=true','-Djosqp.bench.warmupSeconds=3',
                           '-Djosqp.bench.solutionFile='+str(dump),'-cp',str(ROOT/'target'/(mode+'-classes')),
                           'com.quantego.josqp.SolverBenchmark',str(file),str(args.warmup),str(args.samples),'100000']
                print('START',file.name,fork,mode,flush=True)
                run = {'file':file.name,'fork':fork,'mode':mode,'rows':[]}
                start = time.monotonic()
                try:
                    result = subprocess.run(command,capture_output=True,text=True,timeout=args.timeout)
                    stdout = result.stdout
                    run.update(exit_code=result.returncode,stderr=result.stderr[-6000:])
                    if result.returncode: run['failure']='process_error'
                except subprocess.TimeoutExpired as e:
                    stdout = e.stdout or b''
                    run.update(failure='external_timeout',stderr=e.stderr or b'')
                    if isinstance(stdout,bytes): stdout=stdout.decode(errors='replace')
                    if isinstance(run['stderr'],bytes): run['stderr']=run['stderr'].decode(errors='replace')
                run['elapsed_seconds']=time.monotonic()-start
                run['rows']=list(csv.DictReader(io.StringIO(stdout)))
                if run['rows'] and dump.exists(): run['solution_sha256']=hashlib.sha256(dump.read_bytes()).hexdigest()
                data['runs'].append(run)
                output.write_text(json.dumps(data,indent=2)+'\n')
                print('END',file.name,fork,mode,run.get('failure',[(r['solve_ms'],r['allocated_bytes'],r['rho_updates']) for r in run['rows']]),flush=True)
    if any(r.get('failure') for r in data['runs']): raise SystemExit(1)
