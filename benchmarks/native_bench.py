#!/usr/bin/env python3
"""Run native OSQP (and optionally PIQP) over the whole selected case list.

The case list is the same set of historically flagged entries used by
benchmarks/stability.py: every '*' line in the two report tables, minus the
resource cases and missing inputs. Every case is exported once with
ParsedProblemExport so the native solver sees exactly the matrices JOSQP reads.

Example:
  ../.tools/native-env/bin/python benchmarks/native_bench.py \
      --solvers osqp --repeats 3 --output benchmarks/results/upstream-native-bench.json
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from native_compare import read_export, solve

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
SIBLING = WORKSPACE / 'josqp_benchmarks'
EXTRACTED = SIBLING / 'extracted'
JAVA = Path(os.environ.get('JOSQP_JAVA_HOME', os.environ.get('JAVA_HOME',
             '/opt/homebrew/opt/openjdk@26/libexec/openjdk.jdk/Contents/Home')))
BIN_DIR = ROOT / 'target' / 'native-bins'
CLASSES = ROOT / 'target' / 'vector-classes'
EXCLUDED = {'qplib_8500', 'qplib_8547', 'qplib_9008'}
FIELDS = ('status', 'finite', 'setup_ms', 'solve_ms', 'objective',
          'primal_violation', 'stationarity', 'gap_expression')


def selected_cases():
    names = []
    for report in ('report_MarosMeszaros.txt', 'report_QPLIB.txt'):
        for line in (SIBLING / 'Results' / report).read_text().splitlines():
            if line.startswith('*') and '|' in line:
                name = line.split('|')[1].strip()
                if name not in EXCLUDED and name not in names:
                    names.append(name)
    return names


def export_case(name):
    source = EXTRACTED / (name + '.qps')
    target = BIN_DIR / (name + '.native.bin')
    if not source.exists():
        return None, None
    command = [str(JAVA / 'bin/java'), '--add-modules', 'jdk.incubator.vector',
               '-cp', str(CLASSES), 'com.quantego.josqp.ParsedProblemExport',
               str(source), str(target)]
    subprocess.run(command, check=True, capture_output=True)
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    return target, digest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', nargs='+')
    parser.add_argument('--solvers', nargs='+', default=['osqp'])
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--eps', type=float, default=1e-6)
    parser.add_argument('--time-limit', type=float, default=15)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    def case_name(value):
        base = Path(value).name
        return base[:-4] if base.endswith('.qps') else base

    names = [case_name(f) for f in args.files] if args.files else selected_cases()
    BIN_DIR.mkdir(parents=True, exist_ok=True)
    result = {'arguments': vars(args), 'java': str(JAVA), 'cases': [],
              'missing_inputs': [], 'export_failures': []}

    for index, name in enumerate(names):
        binary, digest = export_case(name)
        if binary is None:
            result['missing_inputs'].append(name)
            print('MISSING', name, flush=True)
            continue
        data = read_export(binary)
        case = {'name': name, 'file': str(binary.relative_to(ROOT)), 'n': data['n'],
                'm': data['m'], 'export_sha256': digest, 'solvers': {}}
        for solver in args.solvers:
            entry = {'runs': []}
            for _ in range(args.repeats):
                try:
                    entry['runs'].append({k: solve(data, solver, args.eps,
                                                   time_limit=args.time_limit).get(k) for k in FIELDS})
                except Exception as exc:
                    entry['runs'].append({'failure': type(exc).__name__, 'detail': str(exc)})
            measured = [r for r in entry['runs'] if 'solve_ms' in r]
            if measured:
                entry['median_solve_ms'] = statistics.median(r['solve_ms'] for r in measured)
                entry['median_setup_ms'] = statistics.median(r['setup_ms'] for r in measured)
                entry['statuses'] = sorted({r['status'] for r in measured})
                entry['objectives'] = [r['objective'] for r in measured]
            case['solvers'][solver] = entry
        result['cases'].append(case)
        print('%3d/%d %-22s n=%-7d m=%-7d %s' % (
            index + 1, len(names), name, data['n'], data['m'],
            {s: case['solvers'][s].get('median_solve_ms') for s in args.solvers}), flush=True)

    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
