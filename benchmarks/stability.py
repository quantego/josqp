#!/usr/bin/env python3
"""Audit historically flagged cases, preserving errors and checks in original units."""
import argparse,csv,hashlib,io,json,subprocess,time,tarfile,shutil
from pathlib import Path
from run import ROOT,JAVA,build
BASE='d6eea4c'
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('phase',choices=['before','after','tune'])
 p.add_argument('--files',nargs='+');p.add_argument('--seconds',type=float,default=5)
 p.add_argument('--timeout',type=float,default=35);p.add_argument('--profile',choices=['default','strict','normalized'],default='default')
 p.add_argument('--forks',type=int,default=1);p.add_argument('--output')
 p.add_argument('--eps',type=float);p.add_argument('--eps-rel',type=float);p.add_argument('--max-iter',type=int,default=100000)
 p.add_argument('--no-polish',action='store_true');p.add_argument('--no-gap',action='store_true')
 p.add_argument('--revision');p.add_argument('--linear-audit',action='store_true')
 p.add_argument('--refine',type=int,default=0);p.add_argument('--kkt-rho',action='store_true');p.add_argument('--warmup',action='store_true')
 p.add_argument('--default-rho-interval',action='store_true',help='leave the solver default instead of forcing 100');a=p.parse_args()
 sibling=ROOT.parent/'josqp_benchmarks';out=sibling/'extracted';out.mkdir(exist_ok=True)
 if a.files:names=a.files
 else:
  names=[]
  for report in ['report_MarosMeszaros.txt','report_QPLIB.txt']:
   for line in (sibling/'Results'/report).read_text().splitlines():
    if line.startswith('*') and '|' in line:
     name=line.split('|')[1].strip()+'.qps'
     if name not in ['qplib_8500.qps','qplib_8547.qps','qplib_9008.qps']:names.append(name)
 for name in names:
  dst=out/name;archive=sibling/'QP_Problems'/(Path(name).stem+'.zip')
  if dst.exists() or not archive.exists():continue
  with tarfile.open(archive,'r:*') as tar:
   for m in tar.getmembers():
    if m.isfile() and Path(m.name).name==name:
     with tar.extractfile(m) as src,dst.open('wb') as f:shutil.copyfileobj(src,f)
 revision=a.revision or (BASE if a.phase=='before' else None)
 base=build('vector',revision=revision)
 source_root=ROOT/'target'/('vector-'+revision+'-src') if revision else ROOT
 source_hashes={str(f.relative_to(source_root)):hashlib.sha256(f.read_bytes()).hexdigest() for folder in ['src/main/java','src/vector/java'] for f in sorted((source_root/folder).rglob('*.java'))}
 source_hashes.update({str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((ROOT/'benchmarks/java').rglob('*.java'))})
 result={'phase':a.phase,'arguments':vars(a),'baseline_commit':BASE,'solver_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
  'measured_solver_commit':subprocess.check_output(['git','rev-parse',revision or 'HEAD'],cwd=ROOT,text=True).strip(),
  'source_sha256':source_hashes,
  'java':subprocess.run([str(JAVA/'bin/java'),'-version'],capture_output=True,text=True).stderr,'heap':'2g','runs':[]}
 output=ROOT/(a.output or ('benchmarks/results/stability-'+a.phase+'-'+a.profile+'.json'))
 for name in names:
  file=out/name
  for fork in range(a.forks):
   r={'file':name,'fork':fork,'profile':a.profile,'rows':[]};print('START',name,a.profile,fork,flush=True)
   if not file.exists():r['failure']='missing_input'
   else:
    r['data_sha256']=hashlib.sha256(file.read_bytes()).hexdigest()
    props=['-Daudit.seconds='+str(a.seconds),'-Daudit.maxIter='+str(a.max_iter)]
    if a.profile!='default':props+=['-Daudit.eps=0.000001','-Daudit.infeasEps=0.00000001','-Daudit.polish=true']
    if a.eps is not None:props+=['-Daudit.eps='+str(a.eps)]
    if a.eps_rel is not None:props+=['-Daudit.epsRel='+str(a.eps_rel)]
    if a.no_polish:props+=['-Daudit.polish=false']
    if a.no_gap:props+=['-Daudit.gap=false']
    if a.linear_audit:props+=['-Daudit.linear=true']
    if a.refine:props+=['-Daudit.refine='+str(a.refine)]
    if a.kkt_rho:props+=['-Daudit.kktRho=true']
    if a.warmup:props+=['-Daudit.warmup=true']
    if a.default_rho_interval:props+=['-Daudit.rhoInterval=-1']
    if a.profile=='normalized':
     # Explicit equivalent objective rescaling; divide duals and cost back in audit.
     props+=['-Daudit.costScale=10000']
    cmd=[str(JAVA/'bin/java'),'-Xmx2g','--add-modules','jdk.incubator.vector','-Djosqp.vector=true']+props+['-cp',str(base),'com.quantego.josqp.StabilityBenchmark',str(file)]
    start=time.monotonic()
    try:
     x=subprocess.run(cmd,capture_output=True,text=True,timeout=a.timeout)
     r.update(exit_code=x.returncode,stderr=x.stderr[-5000:]);stdout=x.stdout
     if x.returncode:r['failure']='process_error'
    except subprocess.TimeoutExpired as e:
     r['failure']='external_timeout';stdout=e.stdout or b'';r['stderr']=(e.stderr or b'').decode(errors='replace') if isinstance(e.stderr,bytes) else e.stderr
     if isinstance(stdout,bytes):stdout=stdout.decode(errors='replace')
    r['elapsed_seconds']=time.monotonic()-start
    r['rows']=list(csv.DictReader(io.StringIO(stdout)))
   result['runs'].append(r);output.write_text(json.dumps(result,indent=2)+'\n')
   print('END',name,r.get('failure',[(x['status'],x['objective'],x['primal_violation'],x['stationarity'],x['solve_ms']) for x in r['rows']]),flush=True)
