"""Offline entry point. Run `python -m research.mortality.workbench --help`."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from .runtime import BUNDLE_SHA256, DISCLAIMER, MAX_BODY, MODELS, example, load_bundle, run, schema, strict_json
from .report import write_package


def doctor():
    bundle=load_bundle()
    from .server import check_assets
    check_assets()
    return {'application_version':'0.7.0','runtime_ready':True,
        'static_assets_present':True,'standard_library_only':True,'frozen_model_integrity':'verified',
        'model_sha256':BUNDLE_SHA256,'available_models':list(bundle['models']),
        'historical_US_all_cause_replay_ready':True,'raw_input_imputation':False,
        'clinical_use_ready':False,'Russia_calibrated':False,'Germany_calibrated':False,
        'contemporary_US_calibrated':False,'cause_specific_models_ready':False,
        'twenty_year_risk_ready':False,'upstream_conference_snapshot_imported':False,
        'disclaimer':DISCLAIMER}


def read_payload(path):
    path=Path(path)
    if not path.is_file() or path.stat().st_size>MAX_BODY:raise ValueError('Файл не существует или превышает 1 MiB')
    return strict_json(path.read_bytes())


def batch(path,out):
    """JSONL keeps unit/assay metadata on every row. Invalid rows remain explicit."""
    path=Path(path)
    if not path.is_file() or path.stat().st_size>20*MAX_BODY:raise ValueError('Пакет не существует или слишком велик')
    results=[]
    with path.open('rb') as f:
        for lineno,line in enumerate(f,1):
            if lineno>2000:raise ValueError('Не более 2000 строк за один пакет')
            if not line.strip():raise ValueError('Пустая строка JSONL: '+str(lineno))
            try:r=run(strict_json(line))
            except (ValueError,TypeError,RecursionError):
                r={'status':'blocked','clinical_use_ready':False,'probabilities':None,
                   'errors':[{'field':'line','code':'invalid_JSON_or_integrity_failure'}],
                   'disclaimer':DISCLAIMER}
            results.append({'row':lineno, 'result':r})
    if not results:raise ValueError('Пустой пакет')
    # All integrity errors must fail closed, not merely one row.
    load_bundle()
    summary={'status':'batch_complete_with_blockers' if any(r['result']['status']=='blocked' for r in results) else 'batch_complete_research_only',
       'rows':len(results),'calculated':sum(r['result']['status']=='calculated_research_only' for r in results),
       'blocked':sum(r['result']['status']=='blocked' for r in results),'clinical_use_ready':False,
       'records':results,'disclaimer':DISCLAIMER}
    out=Path(out)
    if out.exists() or out.is_symlink():raise ValueError('Выберите новый каталог для пакета')
    out.mkdir(parents=True)
    try:(out/'batch_results.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    except Exception:
        import shutil;shutil.rmtree(out);raise
    return summary


def main(argv=None):
    p=argparse.ArgumentParser(description='Mortality Workbench 0.7 — локальный исследовательский расчёт, не клинический сервис')
    sub=p.add_subparsers(dest='cmd',required=True)
    sub.add_parser('doctor',help='Проверить контрольную сумму модели и готовность компонентов')
    sub.add_parser('schema',help='Вывести единицы, обязательные поля и ограничения')
    demo=sub.add_parser('demo',help='Синтетический пример → HTML/JSON');demo.add_argument('--out',type=Path,required=True)
    demo.add_argument('--model',choices=list(MODELS),default='M1_routine')
    pred=sub.add_parser('predict',help='Проверить JSON с измерениями → HTML/JSON');pred.add_argument('input',type=Path);pred.add_argument('--out',type=Path,required=True)
    ba=sub.add_parser('batch',help='Пакет JSONL, одна исследовательская запись на строку');ba.add_argument('input',type=Path);ba.add_argument('--out',type=Path,required=True)
    ex=sub.add_parser('example',help='Вывести пример JSON без записи данных');ex.add_argument('--model',choices=list(MODELS),default='M1_routine')
    se=sub.add_parser('serve',help='Локальный браузерный интерфейс на 127.0.0.1');se.add_argument('--port',type=int,default=8765)
    st=sub.add_parser('studio',help='Единый интерфейс моделей 0.16');st.add_argument('--port',type=int,default=8766)
    ev=sub.add_parser('evidence',help='Пересобрать реестр 150 прежних оценок (без нового поиска)');ev.add_argument('--out',type=Path,required=True)
    args=p.parse_args(argv)
    try:
        if args.cmd=='doctor':print(json.dumps(doctor(),ensure_ascii=False,indent=2));return 0
        if args.cmd=='schema':print(json.dumps(schema(),ensure_ascii=False,indent=2));return 0
        if args.cmd=='example':print(json.dumps(example(args.model),ensure_ascii=False,indent=2));return 0
        if args.cmd=='studio':
            from ..studio16.server import serve
            serve(args.port);return 0
        if args.cmd=='serve':
            from .server import serve
            serve(args.port);return 0
        if args.cmd=='evidence':
            if args.out.exists():raise ValueError('Выберите новый каталог для реестра')
            from ..build_evidence_exchange import build
            result=build(Path(__file__).resolve().parents[1],args.out)
            print(json.dumps(result,ensure_ascii=False,indent=2));return 0
        if args.cmd=='batch':
            result=batch(args.input,args.out)
            print(json.dumps({k:v for k,v in result.items() if k!='records'},ensure_ascii=False,indent=2))
            return 2 if result['blocked'] else 0
        payload=example(args.model) if args.cmd=='demo' else read_payload(args.input)
        result=run(payload);write_package(result,args.out)
        if args.cmd=='demo':(args.out/'synthetic_input.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print(json.dumps({'status':result['status'],'report':str(args.out/'report.html'),'clinical_use_ready':False},ensure_ascii=False))
        return 2 if result['status']=='blocked' else 0
    except (ValueError,OSError,KeyError,TypeError,RecursionError) as exc:
        print(json.dumps({'status':'failed','error':str(exc),'clinical_use_ready':False},ensure_ascii=False),file=sys.stderr)
        return 2

if __name__=='__main__':raise SystemExit(main())
