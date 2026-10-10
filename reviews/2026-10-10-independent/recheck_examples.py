"""Run the three verbatim Appendix G scripts unchanged in a scratch directory.
No textbook or release files are written. Requires numpy/scipy/matplotlib.
The wrapper records solver state independently of the book's printed interpretation.
"""
from pathlib import Path
import argparse, hashlib, json, os, re, subprocess, sys

ROOT=Path(__file__).resolve().parents[2]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'tmp/review-20261010/examples')
    out=ap.parse_args().out.resolve();out.mkdir(parents=True,exist_ok=True)
    text=(ROOT/'src/plasma_physics_textbook_v2.tex').read_text(encoding='utf-8')
    blocks=re.findall(r'\\begin\{Verbatim\}[^\n]*\n(.*?)\\end\{Verbatim\}',text,re.S)
    if len(blocks)!=3:raise ValueError('Expected 3 Appendix G scripts; inspect changed source first')
    wrapper='''import json, runpy, sys
import numpy as np
ns=runpy.run_path(sys.argv[1],run_name='__main__')
obj=ns.get('sol')
if obj is not None and hasattr(obj,'rms_residuals'):
    result=dict(status=int(obj.status),success=bool(obj.success),message=obj.message,
                mesh_nodes=len(obj.x),max_rms_residual=float(np.max(obj.rms_residuals)))
    print('REVIEW_SOLVER_STATE='+json.dumps(result))
'''
    wrap=out/'inspect_result.py';wrap.write_text(wrapper,encoding='utf-8')
    env=os.environ.copy();env.update(MPLBACKEND='Agg',PYTHONIOENCODING='utf-8')
    rows=[]
    for i,block in enumerate(blocks,1):
        p=out/f'g{i}.py';p.write_text(block,encoding='utf-8')
        r=subprocess.run([sys.executable,str(wrap),str(p)],cwd=out,env=env,capture_output=True,
                         text=True,encoding='utf-8',timeout=120)
        row=dict(example=i,source_sha256=hashlib.sha256(block.encode()).hexdigest(),
                 exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr.replace(str(out),'$OUT'))
        match=re.search(r'^REVIEW_SOLVER_STATE=(.+)$',r.stdout,re.M)
        if match:row['solver_state']=json.loads(match.group(1))
        rows.append(row)
    (out/'example-runs.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    for row in rows:print(row['example'],row['exit_code'],row.get('solver_state','no BVP result'))
    if any(row['exit_code'] for row in rows):raise SystemExit(1)
    print('Script process completion only; inspect solver_state.success before accepting numerical results.')
if __name__=='__main__':main()
