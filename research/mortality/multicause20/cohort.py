"""Documented fixed-width NHEFS1992/NHANESI inputs, no participant exports.

Only baseline variables predict subsequent outcomes. Native ICD9 categories
remain distinct from the planned ICD10 taxonomy; no national share imputation.
"""
from pathlib import Path
from datetime import date
import hashlib,json,re,calendar
import numpy as np
import pandas as pd

PLAN='858c7e5a7a13799a093d64eb1295b16af2a65a86'
TARGETS=[
 ('ihd','Ишемическая болезнь сердца','410-414','I20-I25'),
 ('heart_failure','Сердечная недостаточность / кардиомиопатии','425,428','I42,I50'),
 ('stroke','Цереброваскулярные заболевания','430-438','I60-I69'),
 ('hypertension','Гипертензивные болезни','401-405','I10-I15'),
 ('atrial_fibrillation','Фибрилляция / трепетание предсердий','4273','I48'),
 ('aortic_aneurysm','Аневризма / диссекция аорты','441','I71'),
 ('chronic_respiratory','Хронические болезни нижних дыхательных путей','490-496','J40-J47'),
 ('dementia','Деменция / болезнь Альцгеймера','290,3310,3312','F01,F03,G30'),
 ('parkinson','Болезнь Паркинсона','332','G20'),
 ('diabetes','Сахарный диабет','250','E10-E14'),
 ('renal','Нефрит / нефроз / почечная недостаточность','580-589','N00-N07,N17-N19,N25-N27'),
 ('chronic_liver','Хроническая болезнь печени / цирроз','571','K70,K73-K74'),
 ('lung_cancer','Рак трахеи, бронхов и лёгкого','162','C33-C34'),
 ('colorectal_cancer','Колоректальный рак','153-154','C18-C21'),
 ('pancreatic_cancer','Рак поджелудочной железы','157','C25'),
 ('breast_cancer','Рак молочной железы','174-175','C50'),
 ('prostate_cancer','Рак простаты','185','C61'),
 ('liver_cancer','Рак печени / внутрипечёночных желчных протоков','155','C22'),
 ('gastric_cancer','Рак желудка','151','C16'),
 ('hematologic_cancer','Гематологические злокачественные опухоли','200-208','C81-C96')]
IDS=[x[0] for x in TARGETS]

def classify(code: str)->str:
    c=str(code).strip().replace('.','')
    if c in {'','0000','000'}:return 'unknown_cause'
    if not re.fullmatch(r'\d{3,4}',c):raise ValueError('ICD9 cause must contain 3/4 digits, implicit E for external')
    n=int(c[:3]);hits=[]
    for key,_,spec,_ in TARGETS:
        for part in spec.split(','):
            if '-' in part:
                lo,hi=map(int,part.split('-'));found=lo<=n<=hi
            else:found=c.startswith(part)
            if found:hits.append(key);break
    if len(hits)>1:raise ValueError('Overlapping target taxonomy')
    if hits:return hits[0]
    if 800<=n<=999:return 'external'
    if 1<=n<=139 or 480<=n<=487:return 'infection'
    # No assertion that this residual excludes all infections elsewhere in ICD9.
    return 'other_or_unresolved'

VITAL={'id':(12,16),'vital':(17,17),'last_m':(18,19),'last_d':(20,21),'last_y':(22,23),
 'exam_m':(33,34),'exam_d':(35,36),'exam_y':(37,38),'age':(39,40),'sex':(46,46),
 'cert':(52,52),'icd9':(53,56),'death_m':(61,62),'death_d':(63,64),'death_y':(65,66),
 'subsample':(219,219),'weight_first35':(220,225),'weight_detail':(232,237),
 'weight':(238,243),'weight_aug':(244,249),'stratum':(259,260),'psu':(261,261)}
LAB={'id':(1,5),'albumin_raw':(232,235),'albumin_flag':(236,236),
 'cholesterol_raw':(237,240),'cholesterol_flag':(241,241),'hemoglobin_raw':(247,250),'hemoglobin_flag':(251,251),
 'bilirubin_raw':(451,454),'ast_raw':(455,458),'alp_raw':(459,462),'urate_raw':(463,465),
 'bun_raw':(472,474),'creatinine_raw':(475,477),'wbc_raw':(529,531)}
BODY={'id':(1,5),'weight_kg_raw':(260,264),'weight_flag':(265,265),'height_raw':(266,269),'height_flag':(273,273)}
BP={'id':(1,5),'sbp':(228,230)}
SMOKE={'id':(1,5),'ever100':(378,378),'current_smoker':(379,379)}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read_fwf(path,fields):
    """Two readers must agree for every selected character field before conversion."""
    raw=Path(path).read_text(encoding='ascii').splitlines()
    width=max(b for a,b in fields.values())
    if any(len(s)<width for s in raw):raise ValueError('Truncated fixed-width record')
    manual=pd.DataFrame({k:[s[a-1:b].strip() for s in raw] for k,(a,b) in fields.items()})
    other=pd.read_fwf(path,colspecs=[(a-1,b) for a,b in fields.values()],names=list(fields),
                      dtype=str,keep_default_na=False,na_filter=False)
    other=other.fillna('').apply(lambda s:s.str.strip())
    if not manual.equals(other):raise ValueError('Independent fixed-width readers disagree')
    if manual.id.duplicated().any() or not manual.id.str.fullmatch(r'\d+').all():raise ValueError('Invalid participant IDs')
    return manual,manual.size

def numeric(s):return pd.to_numeric(s.replace('',np.nan),errors='raise')

def measurement(frame,name,divisor,missing=(777,888,999,7777,8888,9999)):
    v=numeric(frame[name+'_raw']);return (v/divisor).where(~v.isin(missing)&v.gt(0))

def parse_date(y,m,d):
    if not 71<=y<=93 or not 1<=m<=12:return None,False
    day=d if 1<=d<=calendar.monthrange(1900+int(y),int(m))[1] else 15
    return date(1900+int(y),int(m),int(day)),day!=d

def load(source: Path):
    source=Path(source);manifest=json.loads((source/'manifest.json').read_text())
    if manifest['analysis_plan_commit']!=PLAN:raise ValueError('Different input acquisition plan')
    for f in manifest['files']:
        if Path(f['file']).name!=f['file'] or sha(source/f['file'])!=f['sha256']:raise ValueError('Source integrity failure')
    d,cells=read_fwf(source/'n92vitl.txt',VITAL);audit={'vital_n':len(d),'source_hashes_verified':len(manifest['files'])}
    for name,fields in [('DU4800.txt',LAB),('DU4111.txt',BODY),('DU4233.txt',BP),('DU4091.txt',SMOKE)]:
        f,n=read_fwf(source/name,fields);cells+=n;d=d.merge(f,on='id',how='left',validate='1:1')
    # The vital file is the primary endpoint source; certificate file independently crosschecks it.
    cert,n=read_fwf(source/'N92mort.txt',{'id':(3,7),'icd_certificate':(142,145)});cells+=n
    check=d[['id','vital','icd9','cert']].merge(cert,on='id',validate='1:1',how='left')
    known=check.cert.eq('1');mismatch=known&check.icd9.ne(check.icd_certificate)
    if mismatch.any():raise ValueError('Underlying causes differ between source files')
    audit.update(reader_cells_verified=cells,certificate_causes_compared=int(known.sum()),certificate_mismatches=0)
    for col in list(VITAL)+list(BODY)+['sbp','ever100','current_smoker']:
        if col not in {'id','icd9'}:d[col]=numeric(d[col])
    if not d.vital.isin([1,3,4,5,6]).all():raise ValueError('Unrecognized vital status')
    d['dead']=d.vital.eq(3).astype(int);d['cause']=np.where(d.dead.eq(1),d.icd9.map(classify),'alive')
    if ((d.sex==2)&(d.cause=='prostate_cancer')).any():raise ValueError('Sex-incompatible prostate endpoint requires source review')
    years=[];adjust=0;missing_dates=0
    for r in d.itertuples():
        start,a=parse_date(r.exam_y,r.exam_m,r.exam_d)
        end,b=parse_date(r.death_y,r.death_m,r.death_d) if r.dead else parse_date(r.last_y,r.last_m,r.last_d)
        adjust+=int(a)+int(b)
        if start is None or end is None:years.append(np.nan);missing_dates+=1
        else:years.append((end-start).days/365.25)
    d['time']=years
    for name,scale in [('albumin',10),('cholesterol',1),('hemoglobin',10),('creatinine',10),('bun',10),('ast',10),('alp',10),('urate',10),('bilirubin',100),('wbc',10)]:
        d[name]=measurement(d,name,scale)
    for name in ['albumin','cholesterol','hemoglobin']:
        d[name]=d[name].where(numeric(d[name+'_flag']).eq(0))
    mass=measurement(d,'weight_kg',100,missing=(88888,99999));height=measurement(d,'height',10,missing=(8888,9999))
    d['bmi']=(mass/(height/100)**2).where(d.weight_flag.eq(0)&d.height_flag.eq(0))
    d['sbp']=d.sbp.where(d.sbp.gt(0)&~d.sbp.isin([777,888,999]))
    d['smoking_current']=np.where(d.ever100.eq(2),0,np.where(d.ever100.eq(1)&d.current_smoker.isin([1,2]),d.current_smoker.eq(1).astype(float),np.nan))
    d['smoking_former']=np.where(d.ever100.eq(2),0,np.where(d.ever100.eq(1)&d.current_smoker.isin([1,2]),d.current_smoker.eq(2).astype(float),np.nan))
    d['age_scaled']=(d.age-50)/10;d['age_squared']=d.age_scaled**2
    d['albumin']=d.albumin*10 # g/dL -> g/L
    valid=d.age.between(25,74)&d.sex.isin([1,2])&d.time.gt(0)&d.stratum.notna()&d.psu.notna()
    primary=valid&d.subsample.isin([1,2])&d.weight.gt(0)
    audit.update(missing_date_n=missing_dates,day_substitutions=adjust,positive_followup_valid_n=int(valid.sum()),
       primary_eligible_n=int(primary.sum()),source_imputation_excluded={x:int(numeric(d[x+'_flag']).eq(1).sum()) for x in ['albumin','cholesterol','hemoglobin']},
       available_by_subsample={str(int(k)):{f:int(np.isfinite(g[f]).sum()) for f in ['albumin','cholesterol','hemoglobin','creatinine','bun','ast','alp','urate','smoking_current']} for k,g in d.groupby('subsample')})
    return d.loc[primary].reset_index(drop=True),audit

CLINICAL=['age_scaled','age_squared','bmi','sbp']
CORE=CLINICAL+['albumin','cholesterol','hemoglobin']
HEPATIC=CORE+['ast','alp','urate']
