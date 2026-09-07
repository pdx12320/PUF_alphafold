"""Recompute current analyses from the included construct-level caches."""
import os,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
env=os.environ.copy()
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
for name in ['stats.py','models.py','plot_static.py','finish.py','write_report.py']:
 print('Running',name,flush=True)
 subprocess.run([sys.executable,str(root/'architecture_validation'/name)],cwd=root,env=env,check=True)
