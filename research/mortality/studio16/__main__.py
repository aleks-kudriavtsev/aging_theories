"""Unified local research application, standard library runtime."""
import argparse
import json
import sys
from pathlib import Path
from .registry import REGISTRY, catalog
from .engine import run, compare, example, blocked
from .report import save
from ..workbench.runtime import strict_json, MAX_BODY


def read(path):
    path=Path(path)
    if not path.is_file() or path.stat().st_size>MAX_BODY:raise ValueError('Input missing or above 1 MiB')
    return strict_json(path.read_bytes())


def batch(path,out):
    path=Path(path)
    if not path.is_file() or path.stat().st_size>20*MAX_BODY:raise ValueError('Invalid batch size')
    catalog();rows=[]
    with path.open('rb') as f:
        for n,line in enumerate(f,1):
            if n>2000:raise ValueError('At most 2000 rows')
            try:r=run(strict_json(line))
            except (ValueError,TypeError,KeyError,RecursionError,OverflowError):r=blocked([{'field':'line','code':'invalid_JSON_or_artifact'}])
            rows.append({'row':n,'result':r})
    if not rows:raise ValueError('Empty batch')
    catalog() # A damaged artifact cannot be downgraded to a mere row failure.
    d={'version':'0.16.0','status':'batch_research_only','clinical_use_ready':False,'rows':rows,
       'calculated':sum(r['result']['status']=='calculated_research_only' for r in rows),
       'blocked':sum(r['result']['status']=='blocked' for r in rows)}
    out=Path(out)
    if out.exists() or out.is_symlink():raise ValueError('Output exists')
    with out.open('x',encoding='utf-8') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False)
    return d


def main(argv=None):
    p=argparse.ArgumentParser(description='Mortality Studio 0.16 — единый исследовательский интерфейс')
    sub=p.add_subparsers(dest='cmd',required=True)
    sub.add_parser('catalog');sub.add_parser('doctor')
    se=sub.add_parser('serve');se.add_argument('--port',type=int,default=8766)
    ex=sub.add_parser('example');ex.add_argument('--model',choices=list(REGISTRY),required=True)
    de=sub.add_parser('demo');de.add_argument('--model',choices=list(REGISTRY),required=True);de.add_argument('--out',type=Path,required=True)
    pr=sub.add_parser('predict');pr.add_argument('input',type=Path);pr.add_argument('--out',type=Path,required=True)
    co=sub.add_parser('compare');co.add_argument('input',type=Path);co.add_argument('--models',nargs='+',required=True);co.add_argument('--out',type=Path,required=True)
    ba=sub.add_parser('batch');ba.add_argument('input',type=Path);ba.add_argument('--out',type=Path,required=True)
    a=p.parse_args(argv)
    try:
        if a.cmd=='serve':
            from .server import serve
            serve(a.port);return 0
        if a.cmd=='doctor':
            from .server import doctor
            r=doctor();print(json.dumps(r,ensure_ascii=False,indent=2));return 0 if r['status']=='ready_research_only' else 2
        if a.cmd in ('catalog','example'):
            print(json.dumps(catalog() if a.cmd=='catalog' else example(a.model),ensure_ascii=False,indent=2));return 0
        if a.cmd=='batch':
            r=batch(a.input,a.out);print(json.dumps({k:v for k,v in r.items() if k!='rows'}));return 2 if r['blocked'] else 0
        data=example(a.model) if a.cmd=='demo' else read(a.input)
        r=compare(data,a.models) if a.cmd=='compare' else run(data);save(r,a.out)
        print(json.dumps({'status':r['status'],'report':str(a.out/'report.html'),'clinical_use_ready':False},ensure_ascii=False))
        return 2 if r['status']=='blocked' else 0
    except (ValueError,KeyError,TypeError,RecursionError,ArithmeticError,OSError) as e:
        print(json.dumps({'status':'failed','error':str(e),'clinical_use_ready':False},ensure_ascii=False),file=sys.stderr);return 2

if __name__=='__main__':raise SystemExit(main())
