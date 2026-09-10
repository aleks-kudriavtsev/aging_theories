"""Install the wheel into an isolated environment and execute outside source."""
from pathlib import Path
import argparse,subprocess,sys,tempfile,json,os
p=argparse.ArgumentParser();p.add_argument('--wheel',type=Path,required=True);a=p.parse_args()
wheel=a.wheel.resolve()
with tempfile.TemporaryDirectory() as d:
 root=Path(d);subprocess.run([sys.executable,'-m','venv',str(root/'venv')],check=True)
 exe=root/'venv'/('Scripts/python.exe' if os.name=='nt' else 'bin/python')
 subprocess.run([str(exe),'-m','pip','install','--no-index','--no-deps',str(wheel)],check=True)
 commands=[['-m','research.mortality.studio16','doctor'],['-m','research.mortality.studio16','demo','--model','compact4','--out',str(root/'demo')],['-m','research.mortality.workbench','doctor']]
 for cmd in commands:
  done=subprocess.run([str(exe),*cmd],cwd=root,check=True,capture_output=True,text=True)
  payload=json.loads(done.stdout)
  if cmd[-1]=='doctor' and 'studio16' in cmd[1]:assert payload['status']=='ready_research_only'
  print(done.stdout)
 # No source path can satisfy imports in this probe.
 probe=subprocess.run([str(exe),'-c',"import research.mortality.studio16 as m;print(m.__file__)"],cwd=root,capture_output=True,text=True,check=True)
 assert str(root/'venv') in probe.stdout,probe.stdout
 print('ISOLATED_WHEEL_INSTALL_PASSED')
