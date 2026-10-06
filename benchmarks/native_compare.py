#!/usr/bin/env python3
"""Native comparison on JOSQP's exact exported data. Native adapters require installed wheels.

PIQP API: https://predict-epfl.github.io/piqp/interfaces/python/getting_started
OSQP API: https://osqp.org/docs/interfaces/python.html
Use --inspect without dependencies to check the export. No solver is invoked in that mode.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import struct
import time


def read_export(path):
    data=Path(path).read_bytes(); pos=0
    def take(fmt,count=1):
        nonlocal pos
        values=struct.unpack_from('>'+str(count)+fmt,data,pos)
        pos+=struct.calcsize('>'+str(count)+fmt)
        return values
    magic,n,m=take('i',3); assert magic==0x4a515031
    offset=take('d')[0]; matrices=[]
    for rows in (n,m):
        nz=take('i')[0]; ap=take('i',n+1); ai=take('i',nz); ax=take('d',nz)
        assert ap[0]==0 and ap[-1]==nz and all(a<=b for a,b in zip(ap,ap[1:]))
        assert all(0<=r<rows for r in ai)
        matrices.append((rows,n,ap,ai,ax))
    q,l,u=take('d',n),take('d',m),take('d',m)
    assert pos==len(data), 'trailing export data'
    return dict(n=n,m=m,offset=offset,matrices=matrices,q=q,l=l,u=u,
                export_sha256=hashlib.sha256(data).hexdigest())


def solve(d,name,eps):
    import numpy as np
    from scipy import sparse
    upper,a=[sparse.csc_matrix((np.array(ax),np.array(ai),np.array(ap)),shape=(rows,cols))
             for rows,cols,ap,ai,ax in d['matrices']]
    p=upper+upper.T-sparse.diags(upper.diagonal())
    q=np.array(d['q']);l=np.array(d['l']);u=np.array(d['u'])
    l[l<=-1e30]=-np.inf;u[u>=1e30]=np.inf
    if name=='osqp':
        import osqp
        solver=osqp.OSQP();t=time.perf_counter()
        solver.setup(P=upper,q=q,A=a,l=l,u=u,verbose=False,eps_abs=eps,eps_rel=eps,
                     eps_prim_inf=1e-8,eps_dual_inf=1e-8,polishing=True,
                     adaptive_rho_interval=100,max_iter=100000,time_limit=15,check_dualgap=True)
        setup=time.perf_counter()-t;t=time.perf_counter();r=solver.solve()
        elapsed=time.perf_counter()-t;status=r.info.status;x,y=r.x,r.y
    else:
        import piqp
        solver=piqp.SparseSolver();solver.settings.verbose=False
        solver.settings.eps_abs=eps;solver.settings.eps_rel=eps
        equal=np.isfinite(l)&np.isfinite(u)&(l==u)
        unequal=(np.isfinite(l)|np.isfinite(u))&~equal
        t=time.perf_counter()
        solver.setup(P=p.tocsc(),c=q,A=a[equal,:].tocsc(),b=l[equal],
                     G=a[unequal,:].tocsc(),h_l=l[unequal],h_u=u[unequal])
        setup=time.perf_counter()-t;t=time.perf_counter();status=str(solver.solve())
        elapsed=time.perf_counter()-t;r=solver.result;x=r.x
        y=np.zeros(d['m']);y[equal]=r.y;y[unequal]=r.z_u-r.z_l
    if x is None or y is None or not (np.isfinite(x).all() and np.isfinite(y).all()):
        return dict(solver=name,status=status,finite=False,setup_ms=setup*1000,solve_ms=elapsed*1000)
    ax=a@x;px=p@x;dual=px+q+a.T@y
    support=float(np.dot(y[(y>0)&np.isfinite(u)],u[(y>0)&np.isfinite(u)])+
                  np.dot(y[(y<0)&np.isfinite(l)],l[(y<0)&np.isfinite(l)]))
    return dict(solver=name,version=importlib.metadata.version(name),status=status,finite=True,
                setup_ms=setup*1000,solve_ms=elapsed*1000,
                objective=float(0.5*x@px+q@x+d['offset']),
                primal_violation=float(max(0,np.max(l-ax,initial=0),np.max(ax-u,initial=0))),
                stationarity=float(np.max(np.abs(dual),initial=0)),
                gap_expression=float(x@px+q@x+support))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('files',nargs='+');p.add_argument('--inspect',action='store_true')
    p.add_argument('--eps',type=float,default=1e-6);p.add_argument('--output',required=True);a=p.parse_args()
    result={'arguments':vars(a),'adapters_executed':False,'runs':[]}
    for file in a.files:
        d=read_export(file);row=dict(file=file,n=d['n'],m=d['m'],export_sha256=d['export_sha256'])
        row['matrices_nnz']=[len(v[4]) for v in d['matrices']]
        if not a.inspect:
            row['solvers']=[]
            for name in ('osqp','piqp'):
                try:
                    row['solvers'].append(solve(d,name,a.eps));result['adapters_executed']=True
                except ImportError as e: row['solvers'].append(dict(solver=name,failure='missing_dependency',detail=str(e)))
                except Exception as e: row['solvers'].append(dict(solver=name,failure=type(e).__name__,detail=str(e)))
        result['runs'].append(row)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
