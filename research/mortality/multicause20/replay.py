"""Standard-library replay of the historical granular model; never fits on input."""
import argparse,hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).parent
MEASURES={'albumin':('g/L','serum'),'cholesterol':('mg/dL','serum'),'hemoglobin':('g/dL','whole_blood'),
          'ast':('NHANESI_activity','serum'),'alp':('NHANESI_activity','serum'),'urate':('mg/dL','serum')}

def strict_json(raw):
    if len(raw)>1048576:raise ValueError('Input exceeds1MiB')
    def pairs(rows):
        result={}
        for k,v in rows:
            if k in result:raise ValueError('Duplicate JSON key')
            result[k]=v
        return result
    return json.loads(raw,object_pairs_hook=pairs,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite JSON')))

def finite(v):
    try:return type(v) in (int,float) and math.isfinite(v)
    except OverflowError:return False

def load(study,root=ROOT):
    if study not in ('primary','hepatic'):raise ValueError('Explicit primary/hepatic study required')
    index=strict_json((root/'results/MODEL_INDEX.json').read_bytes());row=index[study]
    path=Path(row['path'])
    if path.is_absolute() or '..' in path.parts:raise ValueError('Invalid model path')
    raw=(root/path).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=row['sha256']:raise ValueError('Model integrity failure')
    doc=strict_json(raw)
    if doc['clinical_use_ready'] is not False or doc['validation_design']!='heldout_survey_strata_internal':raise ValueError('Model scope changed')
    return doc,row['sha256']

def score(features,sex,model,horizons):
    pp=model['preprocessing'];z=[];trace=[]
    for i,key in enumerate(pp['features']):
        x=features[key]
        if not finite(x):raise ValueError('Nonfinite or missing feature')
        v=min(max(x,pp['lower'][i]),pp['upper'][i]);z.append((v-pp['mean'][i])/pp['sd'][i])
        trace.append({'feature':key,'value':x,'used_value':v,'clipped':x!=v})
    names=model['cause_order'];rates=[];sexindex=0 if sex=='male' else 1
    for c in names:
        m=model['causes'][c];lp=math.fsum(a*b for a,b in zip(z,m['coefficients']))
        rates.append([0.]*3 if m['female_structural_zero'] and sex=='female' else [math.exp(lp+a) for a in m['baseline_log_hazards'][sexindex]])
    result=[]
    for h in horizons:
        s=1.;cif=[0.]*len(names)
        for j,(lo,hi) in enumerate(zip([0,5,10],[5,10,20])):
            dt=max(0,min(h,hi)-lo);total=math.fsum(r[j] for r in rates);mass=-math.expm1(-total*dt)
            if total>0:
                for k,r in enumerate(rates):cif[k]+=s*mass*r[j]/total
            s*=math.exp(-total*dt)
        if abs(math.fsum(cif)+s-1)>1e-12:raise ArithmeticError('Probability mass failure')
        result.append({'horizon_years':h,'survival':s,'all_cause':math.fsum(cif),'causes':dict(zip(names,cif))})
    return result,trace

def example(study='primary'):
    names=list(MEASURES) if study=='hepatic' else list(MEASURES)[:3]
    values={'albumin':42.,'cholesterol':200.,'hemoglobin':14.,'ast':25.,'alp':70.,'urate':5.}
    return {'schema_version':'0.20.0','study':study,'panel':'laboratory','country':'US',
       'dataset_context':'synthetic_example','acknowledge_research_scope':True,'horizons_years':[10,15],
       'clinical':{'age_years':55,'sex':'male','bmi':27.,'sbp_mmHg':130.,**({'smoking':'former'} if study=='hepatic' else {})},
       'measurements':{k:{'value':values[k],'unit':MEASURES[k][0],'matrix':MEASURES[k][1],'assay_context':'NHANESI_legacy'} for k in names}}

def run(p,root=ROOT):
    keys={'schema_version','study','panel','country','dataset_context','acknowledge_research_scope','horizons_years','clinical','measurements'}
    if not isinstance(p,dict) or set(p)!=keys:raise ValueError('Exact documented fields required; identifiers not accepted')
    if p['schema_version']!='0.20.0' or p['country']!='US' or p['dataset_context'] not in ('synthetic_example','historical_NHEFS_research') or p['acknowledge_research_scope'] is not True:raise ValueError('Historical research scope not confirmed')
    doc,digest=load(p['study'],root);panel=p['panel']
    if panel not in ('clinical','laboratory'):raise ValueError('Explicit panel required')
    c=p['clinical'];required={'age_years','sex','bmi','sbp_mmHg'}|({'smoking'} if p['study']=='hepatic' else set())
    if not isinstance(c,dict) or set(c)!=required:raise ValueError('Incomplete clinical profile')
    if c['sex'] not in ('male','female') or not finite(c['age_years']) or not 25<=c['age_years']<=74:raise ValueError('Unsupported age or sex')
    if not finite(c['bmi']) or not 10<=c['bmi']<=80 or not finite(c['sbp_mmHg']) or not 50<=c['sbp_mmHg']<=300:raise ValueError('Outside software input bounds')
    h=p['horizons_years']
    if not isinstance(h,list) or not h or any(not finite(v) or v not in [10,15] for v in h) or len(set(h))!=len(h):raise ValueError('Only unique10/15year horizons')
    names=set(MEASURES) if p['study']=='hepatic' else {'albumin','cholesterol','hemoglobin'}
    if not isinstance(p['measurements'],dict) or set(p['measurements'])!=names:raise ValueError('Complete same-domain laboratory inputs required')
    x={'age_scaled':(c['age_years']-50)/10,'age_squared':((c['age_years']-50)/10)**2,'bmi':c['bmi'],'sbp':c['sbp_mmHg']}
    for k,r in p['measurements'].items():
        if not isinstance(r,dict) or set(r)!={'value','unit','matrix','assay_context'}:raise ValueError('Measurement metadata missing')
        if (r['unit'],r['matrix'])!=MEASURES[k] or r['assay_context']!='NHANESI_legacy':raise ValueError('No automatic assay or material transfer')
        if not finite(r['value']) or not 0<r['value']<100000:raise ValueError('Invalid measurement; no imputation')
        x[k]=r['value']
    if p['study']=='hepatic':
        if c['smoking'] not in ('never','former','current'):raise ValueError('Smoking must be observed')
        x.update(smoking_current=float(c['smoking']=='current'),smoking_former=float(c['smoking']=='former'))
    predictions,trace=score(x,c['sex'],doc['models'][panel],h)
    selected=doc['estimated_noninfectious_families']
    for r in predictions:
        r['selected_noninfectious_CIF']=math.fsum(r['causes'][k] for k in selected)
        r['unresolved_noninfectious_probability_bounds']={'lower':r['selected_noninfectious_CIF'],
          'upper':r['all_cause']-r['causes']['infection']-r['causes']['external'],
          'kind':'logical_identification_bounds_not_confidence_interval'}
    return {'status':'calculated_historical_research_only','version':'0.20.0','clinical_use_ready':False,
      'study':p['study'],'panel':panel,'model_sha256':digest,'probabilities':predictions,'feature_trace':trace,
      'selected_noninfectious_families':selected,'target_taxonomy_size':20,
      'pooled_unestimated_targets':['atrial_fibrillation','parkinson','liver_cancer'],
      'individual_confidence_interval':None,'imputation_applied':False,'fit_on_request':False,
      'disclaimer':'Историческая когорта NHANES I / NHEFS, 25–74 года. Проверка на исключённых стратах, не внешняя RU/DE. Методы 1970-х не переносимы на новые анализаторы без проверки. Не для клинических решений.'}

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('input',type=Path,nargs='?');a.add_argument('--demo',action='store_true');a.add_argument('--study',choices=['primary','hepatic'],default='primary');a.add_argument('--out',type=Path,required=True)
    args=a.parse_args()
    if bool(args.input)==args.demo:a.error('Input or --demo required, not both')
    if args.out.exists():a.error('Output already exists')
    try:r=run(example(args.study) if args.demo else strict_json(args.input.read_bytes()))
    except (ValueError,KeyError,TypeError,OverflowError,OSError) as e:r={'status':'blocked','error':str(e),'probabilities':None,'clinical_use_ready':False}
    with args.out.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2,allow_nan=False)
    print(json.dumps({'status':r['status'],'clinical_use_ready':False}))
