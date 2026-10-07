#!/usr/bin/env python3
"""Build without Maven using real, locally built JUnit/Hamcrest dependencies.

The current test suite runs on JUnit 5 (Jupiter); it is driven by the JUnit Platform
console launcher, which is fetched from Maven Central into ``.tools/junit5-jars`` on
first use. The archived revisions (baseline, pre-workspace) still contain JUnit 4
sources and keep the old runner.
"""
import argparse
import pathlib
import subprocess
import os
import shutil
import io
import tarfile
import re
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
JAVA = pathlib.Path(os.environ.get('JOSQP_JAVA_HOME', os.environ.get('JAVA_HOME',
                    '/opt/homebrew/opt/openjdk@26/libexec/openjdk.jdk/Contents/Home')))
DEPS = os.environ.get('JOSQP_TEST_CLASSPATH', str(WORKSPACE / '.tools/classes'))
BASELINE = '36f06131eeb8977fa66b31936449802c9401a7cc'
PRE_WORKSPACE = '9a27a42'
JUNIT5_VERSION = '1.14.4'
HAMCREST_VERSION = '2.2'
JUNIT5_DIR = WORKSPACE / '.tools/junit5-jars'
JUNIT5_JAR = JUNIT5_DIR / f'junit-platform-console-standalone-{JUNIT5_VERSION}.jar'
HAMCREST_JAR = JUNIT5_DIR / f'hamcrest-{HAMCREST_VERSION}.jar'
JUNIT4_MODES = ('baseline', 'baseline-qplib', 'pre-workspace')

def fetch(url, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    print(f'Downloading {url}')
    with urllib.request.urlopen(url, timeout=120) as response, open(destination, 'wb') as target:
        shutil.copyfileobj(response, target)

def junit5_classpath():
    if not JUNIT5_JAR.exists():
        fetch('https://repo1.maven.org/maven2/org/junit/platform/junit-platform-console-standalone/'
              f'{JUNIT5_VERSION}/junit-platform-console-standalone-{JUNIT5_VERSION}.jar', JUNIT5_JAR)
    if not HAMCREST_JAR.exists():
        fetch('https://repo1.maven.org/maven2/org/hamcrest/hamcrest/'
              f'{HAMCREST_VERSION}/hamcrest-{HAMCREST_VERSION}.jar', HAMCREST_JAR)
    return os.pathsep.join((str(JUNIT5_JAR), str(HAMCREST_JAR)))

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

def build(mode, revision=None):
    suffix = mode + ('-' + revision if revision else '')
    base = ROOT / 'target' / (suffix + '-classes')
    if base.exists():
        shutil.rmtree(base)
    if revision or mode in ('baseline', 'baseline-qplib', 'pre-workspace'):
        archived_commit = revision or (PRE_WORKSPACE if mode == 'pre-workspace' else BASELINE)
        archived_paths = ['src/main/java', 'src/vector/java'] if revision or mode == 'pre-workspace' else ['src/main/java']
        if not revision:
            archived_paths += ['src/test/java']
            if mode == 'pre-workspace':
                archived_paths += ['src/vector-test/java']
        archive = subprocess.run(['git', 'archive', archived_commit] + archived_paths, cwd=ROOT,
                                 capture_output=True, check=True).stdout
        source_root = ROOT / 'target' / (suffix + '-src' if revision else ('pre-workspace-src' if mode == 'pre-workspace' else 'baseline-src'))
        if source_root.exists():
            shutil.rmtree(source_root)
        source_root.mkdir(parents=True, exist_ok=True)
        with tarfile.open(fileobj=io.BytesIO(archive)) as source:
            # Restrict extraction to regular files/directories, including on Python < 3.12.
            for member in source.getmembers():
                destination = (source_root / member.name).resolve()
                destination.relative_to(source_root.resolve())
                if member.isdir():
                    destination.mkdir(parents=True, exist_ok=True)
                elif member.isfile():
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with source.extractfile(member) as src, destination.open('wb') as dst:
                        shutil.copyfileobj(src, dst)
                else:
                    raise ValueError('Unsupported archive entry: ' + member.name)
    main = ROOT / ('target/baseline-src/src/main/java' if mode in ('baseline', 'baseline-qplib') else 'src/main/java')
    if mode == 'pre-workspace':
        main = ROOT / 'target/pre-workspace-src/src/main/java'
    if revision:
        main = source_root / 'src/main/java'
    if mode == 'baseline-qplib':
        # Original solver; only the BOUNDS token-position fix is shared so the same QPs can be read.
        parser_file = main / 'com/quantego/josqp/Parser.java'
        original = parser_file.read_text()
        current = (ROOT / 'src/main/java/com/quantego/josqp/Parser.java').read_text()
        begin = current.index('\t\tboolean hasValue =')
        end = current.index('\n\t\tString colName', begin)
        old_begin = original.index('\t\tint placeColName = 3;')
        old_end = original.index('\n\t\tString colName', old_begin)
        parser_file.write_text(original[:old_begin] + current[begin:end] + original[old_end:])
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
    if mode not in ('baseline', 'baseline-qplib', 'portable'):
        vector_sources = ROOT / ('target/pre-workspace-src/src/vector/java' if mode == 'pre-workspace' else 'src/vector/java')
        if revision:
            vector_sources = source_root / 'src/vector/java'
        compile_to(base, sorted(vector_sources.rglob('*.java')), str(base), vector=True)
    benchmark_sources = sorted((ROOT / 'benchmarks/java').rglob('*.java'))
    if mode in ('baseline', 'baseline-qplib', 'portable'):
        benchmark_sources = [p for p in benchmark_sources if p.name == 'SolverBenchmark.java']
    compile_to(base, benchmark_sources, str(base), vector=True)
    return base

def junit4_test_source(name):
    """Letzte Fassung einer Testdatei, die noch JUnit 4 importierte (aus der Historie)."""
    relative = f'src/test/java/com/quantego/josqp/{name}'
    revisions = subprocess.run(['git', 'rev-list', 'HEAD', '--', relative], cwd=ROOT,
                               capture_output=True, text=True, check=True).stdout.split()
    for revision in revisions:
        blob = subprocess.run(['git', 'show', f'{revision}:{relative}'], cwd=ROOT,
                              capture_output=True, check=True).stdout
        if b'org.junit.Test' in blob:
            return blob
    raise SystemExit(f'Keine JUnit-4-Fassung von {name} in der Historie gefunden')


def tests(mode):
    base = build(mode)
    test_root = ROOT
    if mode in ('baseline', 'baseline-qplib'):
        test_root = ROOT / 'target/baseline-src'
        # The archived upstream suite mixed Jupiter assertions/annotations with JUnit 4 and does
        # not compile as it stands. These two files are taken from the latest revision in which
        # they still used JUnit 4, so the comparison keeps running against the current JUnit 4
        # classpath.
        for name in ('PrimalDualInfeasibilityTest.java', 'SolveLinSysTest.java'):
            (test_root / 'src/test/java/com/quantego/josqp' / name).write_bytes(junit4_test_source(name))
    elif mode == 'pre-workspace':
        test_root = ROOT / 'target/pre-workspace-src'
    sources = sorted((test_root / 'src/test/java').rglob('*.java'))
    if mode not in ('baseline', 'baseline-qplib', 'portable'):
        sources += sorted((test_root / 'src/vector-test/java').rglob('*.java'))
    vector = mode not in ('baseline', 'baseline-qplib', 'portable')
    if mode in JUNIT4_MODES:
        cp = os.pathsep.join((str(base), str(DEPS)))
        compile_to(base, sources, cp, vector=vector)
        classes = ['com.quantego.josqp.' + p.stem for p in sources if p.stem.endswith('Test')]
        command = [str(JAVA / 'bin/java'), '-cp', cp]
        if vector:
            command += ['--add-modules', 'jdk.incubator.vector']
        if mode == 'pre-workspace':
            command += ['-Djosqp.vector=true']
        subprocess.run(command + ['org.junit.runner.JUnitCore'] + classes, check=True, cwd=ROOT)
        return
    cp = os.pathsep.join((str(base), junit5_classpath()))
    compile_to(base, sources, cp, vector=vector)
    command = [str(JAVA / 'bin/java')]
    if vector:
        command += ['--add-modules', 'jdk.incubator.vector', '-Djosqp.vector=true']
    command += ['-jar', str(JUNIT5_JAR), 'execute', '--disable-ansi-colors',
                '--class-path', cp, '--scan-class-path', '--include-classname', '.*Test',
                '--details=summary']
    subprocess.run(command, check=True, cwd=ROOT)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['build', 'test'])
    parser.add_argument('mode', choices=['baseline', 'baseline-qplib', 'portable', 'scalar', 'vector', 'before-hotloops', 'pre-workspace'])
    args = parser.parse_args()
    tests(args.mode) if args.action == 'test' else build(args.mode)
