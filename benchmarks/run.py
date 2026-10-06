#!/usr/bin/env python3
"""Build without Maven using real, locally built JUnit/Hamcrest dependencies."""
import argparse
import pathlib
import subprocess
import os
import shutil
import io
import tarfile
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
JAVA = pathlib.Path(os.environ.get('JOSQP_JAVA_HOME', os.environ.get('JAVA_HOME',
                    '/opt/homebrew/opt/openjdk@26/libexec/openjdk.jdk/Contents/Home')))
DEPS = os.environ.get('JOSQP_TEST_CLASSPATH', str(WORKSPACE / '.tools/classes'))
BASELINE = '36f06131eeb8977fa66b31936449802c9401a7cc'

def compile_to(destination, sources, classpath=None, vector=False, release=None):
    if not sources:
        return
    destination.mkdir(parents=True, exist_ok=True)
    command = [str(JAVA / 'bin/javac'), '-d', str(destination)]
    if classpath:
        command += ['-cp', classpath]
    if vector:
        command += ['--add-modules', 'jdk.incubator.vector']
    if release:
        command += ['--release', str(release)]
    command += [str(p) for p in sources]
    subprocess.run(command, check=True)

def build(mode):
    base = ROOT / 'target' / (mode + '-classes')
    if base.exists():
        shutil.rmtree(base)
    if mode == 'baseline':
        archive = subprocess.run(['git', 'archive', BASELINE, 'src/main/java'], cwd=ROOT,
                                 capture_output=True, check=True).stdout
        source_root = ROOT / 'target/baseline-src'
        source_root.mkdir(parents=True, exist_ok=True)
        with tarfile.open(fileobj=io.BytesIO(archive)) as source:
            source.extractall(source_root, filter='data')
    main = ROOT / ('target/baseline-src/src/main/java' if mode == 'baseline' else 'src/main/java')
    if mode == 'before-hotloops':
        main = ROOT / 'target/before-hotloops-src'
        if main.exists():
            shutil.rmtree(main)
        shutil.copytree(ROOT / 'src/main/java', main)
        original = subprocess.run(['git', 'show', BASELINE + ':src/main/java/com/quantego/josqp/QDLDL.java'],
                                  cwd=ROOT, capture_output=True, text=True, check=True).stdout
        def method(text, name):
            match = re.search(r'static void ' + name + r'\([^)]*\)\s*\{', text)
            end, depth = match.end(), 1
            while depth:
                depth += (text[end] == '{') - (text[end] == '}')
                end += 1
            return match.start(), end
        path = main / 'com/quantego/josqp/QDLDL.java'
        current = path.read_text()
        for name in ['Lsolve', 'Ltsolve']:
            start, end = method(current, name)
            old_start, old_end = method(original, name)
            current = current[:start] + original[old_start:old_end] + current[end:]
        path.write_text(current)
    compile_to(base, sorted(main.rglob('*.java')), release=8)
    if mode not in ('baseline', 'portable'):
        compile_to(base, sorted((ROOT / 'src/vector/java').rglob('*.java')), str(base), vector=True)
    benchmark_sources = sorted((ROOT / 'benchmarks/java').rglob('*.java'))
    if mode in ('baseline', 'portable'):
        benchmark_sources = [p for p in benchmark_sources if p.name == 'SolverBenchmark.java']
    compile_to(base, benchmark_sources, str(base), vector=True)
    return base

def tests(mode):
    base = build(mode)
    sources = sorted((ROOT / 'src/test/java').rglob('*.java'))
    if mode not in ('baseline', 'portable'):
        sources += sorted((ROOT / 'src/vector-test/java').rglob('*.java'))
    if mode == 'baseline':
        # New regression tests intentionally require the new APIs. Run the original suite only.
        sources = [p for p in sources if p.name not in ('RegressionTest.java', 'VectorKernelsTest.java')]
    cp = os.pathsep.join((str(base), str(DEPS)))
    compile_to(base, sources, cp, vector=mode not in ('baseline', 'portable'))
    classes = ['com.quantego.josqp.' + p.stem for p in sources if p.stem.endswith('Test')]
    command = [str(JAVA / 'bin/java'), '-cp', cp]
    if mode not in ('baseline', 'portable'):
        command += ['--add-modules', 'jdk.incubator.vector']
    if mode == 'vector':
        command += ['-Djosqp.vector=true']
    subprocess.run(command + ['org.junit.runner.JUnitCore'] + classes, check=True, cwd=ROOT)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['build', 'test'])
    parser.add_argument('mode', choices=['baseline', 'portable', 'scalar', 'vector', 'before-hotloops'])
    args = parser.parse_args()
    tests(args.mode) if args.action == 'test' else build(args.mode)
