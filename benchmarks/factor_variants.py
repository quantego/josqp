#!/usr/bin/env python3
"""Independent JVMs compare replay and SIMD blocks with the preceding solver."""
import argparse,csv,hashlib,io,json,random,subprocess,time
from run import ROOT,JAVA,build

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--files',nargs='+',default=['qplib_8495.qps','qplib_8559.qps','qplib_8785.qps'])
 p.add_argument('--forks',type=int,default=3);p.add_argument('--solver',action='store_true')
 p.add_argument('--output',required=True);a=p.parse_args()
 bases={'baseline':build('vector',revision='f9ace75'),'current':build('vector')}
 modes=['baseline','schedule','blocks','blocks-scalar']
 output=ROOT/a.output;result={'arguments':vars(a),'baseline_commit':'f9ace75',
  'solver_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
  'source_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()
    for folder in ['src/main/java','src/vector/java','benchmarks/java'] for f in sorted((ROOT/folder).rglob('*.java'))},'runs':[]}
 rng=random.Random(20261006)
 for name in a.files:
  file=ROOT.parent/'josqp_benchmarks/extracted'/name
  if not file.exists():file=ROOT/'src/test/resources'/name
  for fork in range(a.forks):
   order=modes.copy();rng.shuffle(order)
   for mode in order:
    base=bases['baseline' if mode=='baseline' else 'current']
    cmd=[str(JAVA/'bin/java'),'-Xms256m','-Xmx2g','--add-modules','jdk.incubator.vector',
         '-Djosqp.vector='+str(mode!='blocks-scalar').lower(),
         '-Djosqp.factorSchedule='+str(mode!='baseline').lower(),
         '-Djosqp.factorBlocks='+str(mode.startswith('blocks')).lower(),'-cp',str(base)]
    dump=ROOT/'target'/('upstream-factor-'+file.stem+'-'+mode+'-'+str(fork)+'.bin')
    if a.solver:
     cmd[1:1]=['-Djosqp.bench.warmupSeconds=2','-Djosqp.bench.timeLimit=30','-Djosqp.bench.solutionFile='+str(dump)]
     cmd+=['com.quantego.josqp.SolverBenchmark',str(file),'1','2','100000']
    else:cmd+=['com.quantego.josqp.FactorReplayBenchmark',str(file)]
    print('START',name,fork,mode,flush=True);start=time.monotonic()
    row={'file':name,'fork':fork,'mode':mode,'data_sha256':hashlib.sha256(file.read_bytes()).hexdigest()}
    try:
     r=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
     row.update(exit_code=r.returncode,stderr=r.stderr[-5000:],rows=list(csv.DictReader(io.StringIO(r.stdout))))
     if r.returncode:row['failure']='process_error'
     if a.solver and dump.exists():row['solution_sha256']=hashlib.sha256(dump.read_bytes()).hexdigest()
    except subprocess.TimeoutExpired:row.update(failure='external_timeout',rows=[])
    row['elapsed_seconds']=time.monotonic()-start;result['runs'].append(row)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print('END',name,fork,mode,row.get('failure') or [r.get('factor_ms',r.get('solve_ms')) for r in row['rows']],flush=True)
