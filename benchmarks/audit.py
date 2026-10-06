#!/usr/bin/env python3
"""Demonstrate the regressions on the exact upstream solver, using real JUnit."""
import os
import re
import subprocess
from run import build, compile_to, ROOT, JAVA, DEPS

base = build('baseline')
source = (ROOT / 'src/test/java/com/quantego/josqp/RegressionTest.java').read_text()
start = source.index('    @Test public void inPlaceSolveMatchesPublicCopyAndDoesNotMutateInput()')
end = source.index('    @Test public void productDoesNotSwallowInvalidInput()', start)
# The new in-place API did not exist upstream. All behavioral regressions are unchanged.
source = source[:start] + source[end:]
path = ROOT / 'target/baseline-audit-src/com/quantego/josqp/RegressionTest.java'
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(source)
cp = os.pathsep.join((str(base), str(DEPS)))
compile_to(base, [path], cp, release=8)
log = ROOT / 'benchmarks/results/baseline-audit.txt'
with log.open('w') as handle:
    proc = subprocess.run([str(JAVA / 'bin/java'), '-cp', cp, 'org.junit.runner.JUnitCore',
                           'com.quantego.josqp.RegressionTest'], stdout=handle, stderr=subprocess.STDOUT, cwd=ROOT)
text = log.read_text()
log.write_text(text.rstrip() + '\n')
match = re.search(r'Tests run: (\d+),\s+Failures: (\d+)', text)
if proc.returncode != 1 or not match or match.group(1) != match.group(2):
    raise AssertionError('Unexpected upstream regression result; inspect ' + str(log))
print('Upstream:', match.group(1), 'tests,', match.group(2), 'reproduced failures')
