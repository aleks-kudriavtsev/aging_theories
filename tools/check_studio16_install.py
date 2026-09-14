"""Install outside source, UTF-8 subprocess contract, canonical Windows paths."""
from pathlib import Path
import argparse,subprocess,sys,tempfile,json,os
p=argparse.ArgumentParser();p.add_argument('--wheel',type=Path,required=True);a=p.parse_args()
wheel=a.wheel.resolve()
clean_env=dict(os.environ);clean_env.pop('PYTHONPATH',None);clean_env['PYTHONUTF8']='1'
with tempfile.TemporaryDirectory() as d:
 root=Path(d).resolve();subprocess.run([sys.executable,'-m','venv',str(root/'venv')],check=True)
 exe=root/'venv'/('Scripts/python.exe' if os.name=='nt' else 'bin/python')
 subprocess.run([str(exe),'-m','pip','install','--no-index','--no-deps',str(wheel)],cwd=root,env=clean_env,check=True)
 commands=[['-m','research.mortality.studio16','doctor'],['-m','research.mortality.studio16','demo','--model','compact4','--out',str(root/'demo')],['-m','research.mortality.workbench','doctor']]
 for cmd in commands:
  done=subprocess.run([str(exe),*cmd],cwd=root,env=clean_env,capture_output=True,text=True,encoding='utf-8')
  if done.returncode:raise RuntimeError('Isolated command failed: '+str(cmd)+'; '+done.stderr)
  payload=json.loads(done.stdout)
  if cmd[-1]=='doctor' and 'studio16' in cmd[1]:assert payload['status']=='ready_research_only'
  if cmd[-1]=='doctor' and 'workbench' in cmd[1]:assert payload['runtime_ready']
  print(json.dumps(payload,ensure_ascii=True))
 # Path containment uses canonical paths, not an 8.3-vs-long-name string test.
 probe=subprocess.run([str(exe),'-I','-c',"import research.mortality.studio16 as m;print(m.__file__)"],cwd=root,env=clean_env,capture_output=True,text=True,encoding='utf-8',check=True)
 assert Path(probe.stdout.strip()).resolve().is_relative_to((root/'venv').resolve()),probe.stdout
 print('ISOLATED_WHEEL_INSTALL_PASSED')
