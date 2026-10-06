#!/usr/bin/env python3
"""Bounded exploration of the largest checked-in cases; never count a timeout as solved."""
import json
import subprocess
from run import JAVA, ROOT

results = []
for file in ['nug08-3rd.mps', 'irish-electricity.mps', 'ex10.mps', 'supportcase10.mps', 'rmine15_lp.mps']:
    command = [str(JAVA / 'bin/java'), '-Xms256m', '-Xmx2g', '--add-modules', 'jdk.incubator.vector',
               '-Djosqp.vector=true', '-Djosqp.bench.timeLimit=20', '-cp', str(ROOT / 'target/vector-classes'),
               'com.quantego.josqp.SolverBenchmark', str(ROOT / 'src/test/resources' / file), '0', '1', '10000']
    print(file, flush=True)
    try:
        proc = subprocess.run(command, capture_output=True, text=True, timeout=60, cwd=ROOT)
        row = dict(file=file, returncode=proc.returncode, stdout=proc.stdout, stderr=proc.stderr)
    except subprocess.TimeoutExpired as e:
        def decoded(value):
            return value.decode(errors='replace') if isinstance(value, bytes) else value
        row = dict(file=file, timeout_seconds=60, stdout=decoded(e.stdout), stderr=decoded(e.stderr))
    results.append(row)
    (ROOT / 'benchmarks/results/large-probes.json').write_text(json.dumps(results, indent=2))
    print(row, flush=True)
