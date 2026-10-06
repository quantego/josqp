#!/usr/bin/env python3
"""Serial QPLIB pilot and three-variant comparison; retain failures and provenance."""
import argparse
import csv
import hashlib
import io
import json
import random
import subprocess
import time
import tarfile
import shutil
from pathlib import Path
from run import ROOT, JAVA, build, BASELINE

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('phase', choices=['pilot', 'compare'])
    p.add_argument('--timeout', type=int, default=60)
    args = p.parse_args()
    sibling = ROOT.parent / 'josqp_benchmarks'
    extracted = sibling / 'extracted'
    extracted.mkdir(exist_ok=True)
    # The files named .zip in this repository are actually gzip-compressed tar archives.
    for archive in sorted((sibling / 'QP_Problems').glob('qplib_*.zip')):
        with tarfile.open(archive, 'r:gz') as tar:
            for member in tar.getmembers():
                name = Path(member.name).name
                if member.isfile() and name.startswith('qplib_') and name.endswith('.qps'):
                    destination = extracted / name
                    if not destination.exists():
                        with tar.extractfile(member) as src, destination.open('wb') as dst:
                            shutil.copyfileobj(src, dst)
    results = ROOT / 'benchmarks/results'
    files = sorted((sibling / 'extracted').glob('qplib_*.qps'))
    modes = ['vector'] if args.phase == 'pilot' else ['baseline-qplib', 'scalar', 'vector']
    if args.phase == 'compare':
        pilot = json.loads((results / 'qplib-pilot.json').read_text())
        measured = {r['file'] for r in pilot['runs'] if r['rows']}
        files = [f for f in files if f.name in measured]
    for mode in modes:
        build(mode)
    data = {'phase': args.phase, 'max_iter': 100000, 'heap': '2g',
            'timeout_seconds_per_jvm': args.timeout, 'baseline_commit': BASELINE,
            'baseline_parser_patch': 'Only BOUNDS variable-position fix shared; original solver otherwise unchanged',
            'solver_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'data_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=sibling, text=True).strip(),
            'java': subprocess.run([str(JAVA/'bin/java'), '-version'], capture_output=True, text=True).stderr,
            'source_sha256': {str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest()
                              for directory in ['src/main/java', 'src/vector/java', 'benchmarks/java']
                              for f in sorted((ROOT / directory).rglob('*.java'))},
            'data_sha256': {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
            'runs': []}
    output = results / ('qplib-' + args.phase + '.json')
    rng = random.Random(20261006)
    forks = 1 if args.phase == 'pilot' else 3
    for f in files:
        for fork in range(forks):
            order = list(modes); rng.shuffle(order)
            for mode in order:
                cmd = [str(JAVA/'bin/java'), '-Xms256m', '-Xmx2g', '--add-modules', 'jdk.incubator.vector',
                       '-Djosqp.bench.warmupSeconds=' + ('0' if args.phase == 'pilot' else '3'),
                       '-cp', str(ROOT/'target'/(mode+'-classes'))]
                if mode == 'vector': cmd += ['-Djosqp.vector=true']
                dump = ROOT / 'target' / ('qplib-' + f.stem + '-' + mode + '-' + str(fork) + '.bin')
                if args.phase == 'compare':
                    cmd += ['-Djosqp.bench.solutionFile=' + str(dump)]
                cmd += ['com.quantego.josqp.SolverBenchmark', str(f),
                        '0' if args.phase == 'pilot' else '1', '1' if args.phase == 'pilot' else '2', '100000']
                r = {'file': f.name, 'fork': fork, 'mode': mode, 'rows': []}
                start = time.monotonic()
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
                    r.update(exit_code=result.returncode, stderr=result.stderr[-6000:])
                    stdout = result.stdout
                    if result.returncode: r['failure'] = 'process_error'
                except subprocess.TimeoutExpired as e:
                    r['failure'] = 'external_timeout'
                    stdout = e.stdout or b''
                    r['stderr'] = e.stderr or b''
                    if isinstance(stdout, bytes): stdout = stdout.decode(errors='replace')
                    if isinstance(r['stderr'], bytes): r['stderr'] = r['stderr'].decode(errors='replace')
                r['elapsed_seconds'] = time.monotonic() - start
                r['rows'] = list(csv.DictReader(io.StringIO(stdout)))
                if args.phase == 'compare' and dump.exists():
                    r['solution_sha256'] = hashlib.sha256(dump.read_bytes()).hexdigest()
                data['runs'].append(r)
                output.write_text(json.dumps(data, indent=2) + '\n')
                print(f.name, fork, mode, r.get('failure', [(x['status'], x['iterations'], x['solve_ms']) for x in r['rows']]), flush=True)
