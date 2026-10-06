#!/usr/bin/env python3
"""Control runs for the infeasibility diagnosis.

Native OSQP and JOSQP are run on exactly the same exported matrices, first with the JOSQP default
profile (eps_abs = eps_rel = 1e-3, eps_prim_inf = eps_dual_inf = 1e-4, no polishing, adaptive rho
interval 50) and then with a sweep over the infeasibility tolerances. Native OSQP uses the compiled
wheel, so the comparison separates "port defect" from "upstream algorithm behaviour".

Example:
  ../.tools/native-env/bin/python benchmarks/infeasibility_control.py \
      --output benchmarks/results/infeasibility-native-control.json
"""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from native_compare import read_export

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
EXTRACTED = WORKSPACE / 'josqp_benchmarks' / 'extracted'
JAVA = Path(os.environ.get('JOSQP_JAVA_HOME', os.environ.get('JAVA_HOME',
             '/opt/homebrew/opt/openjdk@26/libexec/openjdk.jdk/Contents/Home')))
BINS = ROOT / 'target' / 'native-bins'
CLASSES = ROOT / 'target' / 'vector-classes'
CASES = ('hues-mod', 'primalc5', 'primalc8')


def export(name):
    target = BINS / (name + '.native.bin')
    if not target.exists():
        BINS.mkdir(parents=True, exist_ok=True)
        subprocess.run([str(JAVA / 'bin/java'), '--add-modules', 'jdk.incubator.vector',
                        '-cp', str(CLASSES), 'com.quantego.josqp.ParsedProblemExport',
                        str(EXTRACTED / (name + '.qps')), str(target)], check=True)
    return target


def native(name, eps, eps_inf, interval=50, time_limit=30):
    import numpy as np
    from scipy import sparse
    import osqp
    d = read_export(export(name))
    upper, a = [sparse.csc_matrix((np.array(ax), np.array(ai), np.array(ap)), shape=(rows, cols))
                for rows, cols, ap, ai, ax in d['matrices']]
    q = np.array(d['q'])
    l = np.array(d['l'])
    u = np.array(d['u'])
    l[l <= -1e30] = -np.inf
    u[u >= 1e30] = np.inf
    solver = osqp.OSQP()
    solver.setup(P=upper, q=q, A=a, l=l, u=u, verbose=False, eps_abs=eps, eps_rel=eps,
                 eps_prim_inf=eps_inf, eps_dual_inf=eps_inf, polishing=False,
                 adaptive_rho_interval=interval, max_iter=100000, time_limit=time_limit,
                 check_dualgap=True)
    r = solver.solve()
    return {'solver': 'osqp-native', 'version': osqp.__version__, 'status': r.info.status,
            'iterations': int(r.info.iter), 'objective': float(r.info.obj_val),
            'pri_res': float(r.info.prim_res), 'dua_res': float(r.info.dual_res),
            'eps': eps, 'eps_inf': eps_inf, 'rho_interval': interval}


def josqp(name, streak, margin=None, seconds=30):
    props = ['-Daudit.seconds=' + str(seconds), '-Daudit.maxIter=100000', '-Daudit.rhoInterval=-1',
             '-Daudit.infeasStreak=' + str(streak)]
    if margin is not None:
        props.append('-Daudit.infeasMargin=' + str(margin))
    command = [str(JAVA / 'bin/java'), '-Xmx2g', '--add-modules', 'jdk.incubator.vector',
               '-Djosqp.vector=true'] + props + ['-cp', str(CLASSES),
               'com.quantego.josqp.StabilityBenchmark', str(EXTRACTED / (name + '.qps'))]
    out = subprocess.run(command, capture_output=True, text=True, check=True).stdout
    rows = list(csv.DictReader(io.StringIO(out)))
    row = rows[-1]
    return {'solver': 'josqp', 'infeasibility_streak': streak, 'infeasibility_margin': margin,
            'status': row['status'],
            'iterations': int(row['iterations']), 'objective': float(row['objective']),
            'primal_violation': float(row['primal_violation']),
            'stationarity': float(row['stationarity']), 'solve_ms': float(row['solve_ms'])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    parser.add_argument('--sweep', nargs='+', type=float, default=[1e-3, 1e-4, 1e-5, 1e-6])
    args = parser.parse_args()
    result = {'java': str(JAVA), 'cases': {}}
    for name in CASES:
        binary = export(name)
        entry = {'export_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
                 'default_profile': {
                     'native': native(name, 1e-3, 1e-4),
                     'josqp_default': josqp(name, 2),
                     'josqp_upstream': josqp(name, 1, 1.0)},
                 'tolerance_sweep': [native(name, 1e-3, e) for e in args.sweep]}
        result['cases'][name] = entry
        print(name, json.dumps(entry, indent=1)[:1200], flush=True)
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
