#!/usr/bin/env python3
"""Run independent JVM forks sequentially; save all samples and solver outputs."""
import argparse
import csv
import io
import json
import pathlib
import subprocess
import sys
import random
import hashlib
import platform
from run import ROOT, JAVA, build

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', nargs='+', default=['qafiro.qps', 'qbandm.qps', 'qsctap1.qps', 'qap15.mps'])
    parser.add_argument('--modes', nargs='+', default=['baseline', 'scalar', 'vector'],
                        choices=['baseline', 'scalar', 'vector', 'before-hotloops'])
    parser.add_argument('--forks', type=int, default=3)
    parser.add_argument('--warmup', type=int, default=5)
    parser.add_argument('--warmup-seconds', type=float, default=3)
    parser.add_argument('--samples', type=int, default=7)
    parser.add_argument('--max-iter', type=int, default=10000)
    parser.add_argument('--timeout', type=int, default=120)
    parser.add_argument('--output', default='benchmarks/results/solvers.csv')
    args = parser.parse_args()
    for mode in args.modes:
        build(mode)
    output = ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    source_paths = sorted((ROOT / 'src/main/java').rglob('*.java')) + sorted((ROOT / 'src/vector/java').rglob('*.java'))
    output.with_suffix('.metadata.json').write_text(json.dumps({
        'arguments': vars(args), 'platform': platform.platform(), 'architecture': platform.machine(),
        'java': subprocess.run([str(JAVA / 'bin/java'), '-version'], capture_output=True, text=True).stderr,
        'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
    }, indent=2))
    writer = None
    failures = []
    rng = random.Random(9273)
    with output.open('w') as handle:
        for file in args.files:
            for fork in range(args.forks):
                modes = list(args.modes)
                rng.shuffle(modes)
                for mode in modes:
                    command = [str(JAVA / 'bin/java'), '-Xms256m', '-Xmx2g', '--add-modules', 'jdk.incubator.vector',
                               '-Djosqp.bench.warmupSeconds=' + str(args.warmup_seconds),
                               '-cp', str(ROOT / 'target' / (mode + '-classes'))]
                    if mode == 'vector':
                        command += ['-Djosqp.vector=true']
                    command += ['com.quantego.josqp.SolverBenchmark', str(ROOT / 'src/test/resources' / file),
                                str(args.warmup), str(args.samples), str(args.max_iter)]
                    print(file, 'fork', fork, mode, flush=True)
                    try:
                        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=True, timeout=args.timeout)
                    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
                        failures.append({'file': file, 'fork': fork, 'mode': mode, 'error': str(e),
                                         'stderr': str(getattr(e, 'stderr', ''))[-4000:]})
                        print('failed:', type(e).__name__, flush=True)
                        continue
                    for row in csv.DictReader(io.StringIO(result.stdout)):
                        row['file'] = file
                        row.update(mode=mode, fork=fork)
                        if writer is None:
                            writer = csv.DictWriter(handle, fieldnames=list(row), lineterminator='\n')
                            writer.writeheader()
                        writer.writerow(row)
                        handle.flush()
    output.with_suffix('.failures.json').write_text(json.dumps(failures, indent=2))
    if failures:
        sys.exit(1)
