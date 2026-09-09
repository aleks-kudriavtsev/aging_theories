"""Smoke the installed wheel outside the source tree; synthetic records only."""
import argparse,json,os,subprocess,sys,tempfile
from pathlib import Path

def check(wheel):
    wheel=Path(wheel).resolve()
    with tempfile.TemporaryDirectory(prefix='workbench-install-') as directory:
        root=Path(directory);env=root/'venv'
        subprocess.run([sys.executable,'-m','venv',str(env)],check=True)
        py=env/('Scripts/python.exe' if os.name=='nt' else 'bin/python')
        cli=env/('Scripts/mortality-workbench.exe' if os.name=='nt' else 'bin/mortality-workbench')
        subprocess.run([str(py),'-m','pip','install','--no-index','--no-deps',str(wheel)],cwd=root,check=True)
        clean_env=dict(os.environ);clean_env.pop('PYTHONPATH',None);clean_env['PYTHONUTF8']='1'
        def call(*args):
            result=subprocess.run([str(cli),*args],cwd=root,env=clean_env,capture_output=True,text=True,encoding='utf-8',check=True)
            return json.loads(result.stdout)
        doctor=call('doctor');assert doctor['runtime_ready'] and doctor['static_assets_present'] and not doctor['clinical_use_ready']
        demo=call('demo','--out',str(root/'demo'));assert demo['status']=='calculated_research_only'
        report=json.loads((root/'demo/report.json').read_text(encoding='utf-8'));assert len(report['probabilities'])==3
        evidence=call('evidence','--out',str(root/'evidence'));assert evidence['records']==150
        import_check=subprocess.run([str(py),'-I','-c','from research.mortality.workbench.server import readiness; assert readiness()["status"]=="ready_research_only"'],cwd=root,env=clean_env,check=True)
        return {'installed_outside_source':True,'doctor':True,'demo':True,'evidence_records':150,'isolated_import':True,'clinical_use_ready':False}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('wheel',type=Path);a=p.parse_args();print(json.dumps(check(a.wheel),indent=2))
