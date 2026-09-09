"""Loss-aware migration of frozen evidence and NHANES metrics; no model fitting.

python -m research.mortality.build_evidence_exchange --root research/mortality --out /tmp/exchange
An optional --pr115-snapshot accepts the upstream prediction_performance.json.
Source texts are not silently reverified; migration is not a literature review.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
from .evidence_contract import SCHEMA_VERSION, audit_record, digest, project_estimates

BASE_COMMIT = 'fe414dc1957d8c6041dce48e77cf58ee65b30d28'
UPSTREAM_COMMIT = '3c487a8d03b15a49ce1d067e8eba903dab826b10'
ENDPOINTS = {
 'all_cause_mortality':'all_cause_mortality', 'all_cause_death':'all_cause_mortality',
 'all_cause_mortality_in_established_HF':'all_cause_mortality',
 'respiratory_mortality':'respiratory_mortality', 'cardiovascular_death':'cardiovascular_mortality',
 'liver_related_death':'liver_mortality',
 'incident_heart_failure_NOT_death':'incident_heart_failure',
 'CHD_or_stroke':'CHD_or_stroke_events', 'CHD_or_stroke_or_HF':'CHD_or_stroke_or_HF_events',
 'cardiovascular_events_NOT_all_cause_mortality':'cardiovascular_events',
 'incident_CVD_NOT_all_cause_mortality':'cardiovascular_events',
 'kidney_failure_dialysis_or_transplant_NOT_death':'kidney_failure_RRT',
 'lung_cancer_diagnosis':'incident_lung_cancer',
 'colorectal_cancer':'diagnosis_colorectal_cancer',
 'advanced_precancerous_lesion':'diagnosis_advanced_precancerous_lesion',
 'advanced_colorectal_neoplasia':'diagnosis_advanced_colorectal_neoplasia',
 'Alzheimer_pathology':'diagnosis_Alzheimer_pathology',
 # The legacy unspecific string 'mortality' is intentionally unresolved.
}
METRIC_ALIASES = {'C':'C_index','median_external_C':'C_index',
                 'median_calibration_slope':'calibration_slope'}
COUNTRIES = {'USA':'USA', 'Germany':'DEU', 'Korea':'KOR', 'Sweden':'SWE'}
MODEL_STUDIES = {'PREVENT':'PREVENT', 'KFRE2016':'KFRE',
                'MAGGIC_VALIDATION2014':'MAGGIC'}


def source_ref(path: Path, *, locator: str, kind: str) -> dict:
    result = {'kind':kind, 'repository':'aleks-kudriavtsev/aging_theories',
            'commit':BASE_COMMIT, 'path':'research/mortality/'+('nhanes04/' if path.parent.name=='nhanes04' else '')+path.name,
            'sha256':hashlib.sha256(path.read_bytes()).hexdigest(), 'locator':locator,
            'review':'frozen_empirical_output' if kind=='empirical_evaluation' else 'legacy_extraction_not_reverified'}
    if path.name == 'validation_intervals.csv':
        # The frozen companion archive contains this table; baseline fe414dc
        # did not publish it. Never invent a Git URL at that historical commit.
        result.update(commit=None, code_commit=BASE_COMMIT,
          original_archive='mortality_iteration04.zip',
          original_archive_member='mortality_iteration04/research/mortality/nhanes04/validation_intervals.csv',
          original_archive_sha256='324b76135f7e4477dba0a494bfdad58f7d27613b62c93af8435cdc893f2c6ee7')
    return result


def measurement(label: str | None) -> dict:
    if label in {'eGFR','UACR','FIB4'}:
        kind = 'derived_lab_measure'
    elif label in {'BODE_score','PLCOm2012','PLCOm2012_selection','NLST_selection'}:
        kind = 'clinical_model_or_score'
    elif label in {'gait_speed_m_s','gait_speed_slope'}:
        kind = 'functional_measure'
    elif label == 'FEV1':
        kind = 'physiological_measure'
    elif label and ('MELD' in label):
        kind = 'clinical_model_or_score'
    elif label and '+' in label:
        kind = 'joint_component_block'
    elif label in {'APS2','cfDNA_test'}:
        kind = 'diagnostic_algorithm'
    else:
        kind = 'laboratory_analyte' if label else 'not_reported'
    return {'id':label, 'label':label, 'kind':kind, 'matrix':None, 'assay':None, 'units':None,
            'contains_age':label in {'FIB4','eGFR'}}


def literature_records(root: Path) -> list[dict]:
    rows = []
    for filename in ('evidence.json','evidence_02.json'):
        path = root/filename
        doc = json.loads(path.read_text())
        for si, study in enumerate(doc['studies']):
            for ei, raw in enumerate(study['estimates']):
                endpoint_original = raw.get('outcome',study.get('outcome'))
                metric = METRIC_ALIASES.get(raw['metric'],raw['metric'])
                source = source_ref(path,locator=f'/studies/{si}/estimates/{ei}',kind='literature_extraction')
                source.update(pmid=study.get('pmid'),doi=study.get('doi'),url=study.get('url'),
                              quarantined=study.get('quality_gate')=='quarantine_for_model_training')
                interval = None
                limits = raw.get('ci95',raw.get('interval'))
                if limits is not None:
                    itype = '95% CI' if raw.get('ci95') is not None else raw.get('interval_type')
                    itype = {'interquartile_interval':'IQR'}.get(itype,itype)
                    interval = {'lower':limits[0], 'upper':limits[1],
                                'kind':itype or 'unknown_reported_interval',
                                'scope':'as_reported_in_legacy_extraction'}
                horizon = raw.get('horizon_years',study.get('prediction_horizon_years',study.get('horizon_years')))
                context = {'endpoint':ENDPOINTS.get(endpoint_original), 'endpoint_original':endpoint_original,
                    'population':study.get('population'), 'country':COUNTRIES.get(study.get('country')),
                    'sex':raw.get('sex'), 'age_range':study.get('age_range'),
                    'horizon_years':horizon,'horizon_basis':'prediction_horizon' if horizon is not None else None,
                    'cohort_id':raw.get('sample'), 'validation_set':raw.get('sample'),
                    'sex_age_handling_original':study.get('sex_age_handling',study.get('sex_handling')),
                    'followup_original':{k:v for k,v in study.items() if 'followup' in k},
                    'other_horizon_units':{k:v for k,v in study.items() if k in {'prediction_horizon_days','prediction_horizon_months'}}}
                label = raw.get('marker')
                model_id = MODEL_STUDIES.get(study['id'])
                if measurement(label)['kind'] in {'clinical_model_or_score','diagnostic_algorithm'}:
                    model_id = label
                record = {'schema_version':SCHEMA_VERSION,'id':f"literature:{study['id']}:{ei}",
                    'record_type':'association' if metric in {'HR','RR','OR','reported_relative_reduction','reported_relative_increase'}
                                  else 'diagnostic_accuracy' if raw.get('endpoint_class')=='diagnosis'
                                  else 'prediction_performance',
                    'context':context,'measurement':measurement(label),
                    'model':{'id':model_id, 'base_id':None, 'added_components':[label] if label else [],
                             'component_kind':'marker_block' if label and '+' in label else 'not_extracted'},
                    'estimate':{'metric':metric,'value':raw['value'],'interval':interval,
                                'contrast':raw.get('contrast'), 'summary_statistic':'median_across_validations' if raw['metric'].startswith('median_') else None},
                    'source':source,'raw_record':deepcopy(raw),
                    'clinical_use_ready':False,'usable_as_mortality_coefficient':False}
                rows.append(record)
    return rows


def csv_rows(path: Path) -> list[dict]:
    with path.open(newline='',encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def empirical_records(root: Path) -> list[dict]:
    folder = root/'nhanes04'
    specs = {'M0_age_sex':[], 'M0c_clinical_only_posthoc':[],
             'M1_routine':[], 'M2_NTproBNP':['NT-proBNP'], 'M3_troponin':['hs-cTnT'],
             'M4_cystatinC':['cystatin_C']}
    base = csv_rows(folder/'metrics.csv')
    clinical = [r for r in csv_rows(folder/'clinical_sensitivity.csv') if r['model']=='M0c_clinical_only_posthoc']
    pairs = {(r['model'],int(r['horizon_years'])):r for r in base+clinical}
    results=[]
    fields = {'auc':'fixed_horizon_AUC','brier':'Brier_score',
              'calibration_slope':'calibration_slope',
              'calibration_in_the_large_logit':'calibration_in_the_large',
              'observed_expected':'observed_expected_ratio'}
    def make(row,field,value,locator,filename,comparison=None):
        src=source_ref(folder/filename,locator=locator,kind='empirical_evaluation')
        src['path']='research/mortality/nhanes04/'+filename
        src['pmid']=None  # A local calculation is not a published clinical study.
        context={'endpoint':'all_cause_mortality',
          'population':'NHANES examination; age 40-79; eligible protein-weighted analytic subset',
          'country':'USA','sex':'pooled','age_range':[40,79],
          'horizon_years':int(row['horizon_years']),'horizon_basis':'prediction_horizon',
          'cohort_id':'NHANES_2003_2004_protein_subset_n2384', 'validation_set':'temporal_validation',
          'n':int(row['n']),'events':int(row['events'])}
        model={'id':row['model'],'base_id':None,'added_components':specs.get(row['model'],[]),
               'component_kind':'model', 'training_cohort':'NHANES_1999_2002_n4552',
               'analysis_status':'previously_examined_frozen_benchmark',
               'post_hoc':row['model']=='M0c_clinical_only_posthoc'}
        if comparison:model.update(comparison)
        return {'schema_version':SCHEMA_VERSION,'id':f"empirical:{row['model']}:{row['horizon_years']}:{field}:"+(model['base_id'] or 'absolute'),
          'record_type':'prediction_performance','context':context,'measurement':{},'model':model,
          'estimate':{'metric':field,'value':value,'interval':None}, 'source':src,
          'clinical_use_ready':False,'usable_as_mortality_coefficient':False}
    for row in base+clinical:
        fn='clinical_sensitivity.csv' if row['model']=='M0c_clinical_only_posthoc' else 'metrics.csv'
        for column,metric in fields.items():
            results.append(make(row,metric,float(row[column]),f"model={row['model']};horizon={row['horizon_years']};column={column}",fn))
    # All comparisons are old, explicitly frozen; no new validation is claimed.
    contrasts=[('M1_routine','M0_age_sex',['clinical_and_routine_laboratory_block'],'mixed_block',False),
      ('M1_routine','M0c_clinical_only_posthoc',['routine_laboratory_block'],'marker_block',True),
      ('M2_NTproBNP','M1_routine',['NT-proBNP'],'single_marker',False),
      ('M3_troponin','M2_NTproBNP',['hs-cTnT'],'single_marker',False),
      ('M4_cystatinC','M3_troponin',['cystatin_C'],'single_marker',False),
      ('M4_cystatinC','M1_routine',['NT-proBNP','hs-cTnT','cystatin_C'],'marker_block',False)]
    interval_path=folder/'validation_intervals.csv'
    intervals=csv_rows(interval_path) if interval_path.exists() else []
    ci_path=folder/'clinical_sensitivity_intervals.csv'
    clin_intervals=csv_rows(ci_path)
    for model,base_id,components,kind,posthoc in contrasts:
        row, reference=pairs[(model,10)],pairs[(base_id,10)]
        for col,metric in [('auc','delta_AUC'),('brier','delta_Brier')]:
            rec=make(row,metric,float(row[col])-float(reference[col]),f'{model} minus {base_id};horizon=10;column={col}',
                     'metrics.csv',{'base_id':base_id,'added_components':components,'component_kind':kind,'post_hoc':posthoc})
            if posthoc:
                ref=next(r for r in clin_intervals if r['metric']==f'delta_{col}_routine_minus_clinical')
                ipath=ci_path
                rec['source']['secondary_sources']=[source_ref(folder/'clinical_sensitivity.csv',locator='M0c_clinical_only_posthoc;horizon=10',kind='empirical_evaluation')]
            else:
                name=f'delta_{col}_routine' if (model,base_id)==('M4_cystatinC','M1_routine') else f'delta_{col}_previous'
                ref=next((r for r in intervals if r['model']==model and r['metric']==name),None)
                ipath=interval_path
                # The published repo includes the M4-vs-M1 interval in the clinical CSV.
                if ref is None and (model,base_id)==('M4_cystatinC','M1_routine') and col=='auc':
                    ref=next(r for r in clin_intervals if r['metric']=='delta_auc_expanded_minus_routine');ipath=ci_path
            if ref:
                rec['estimate']['interval']={'kind':'conditional_validation_percentiles',
                  'lower':float(ref['p025']),'upper':float(ref['p975']),'replicates':int(ref['replicates']),
                  'scope':'fixed_models_validation_PSU_resampling; excludes_training_and_model_selection',
                  'source':source_ref(ipath,locator=ref['metric'],kind='empirical_evaluation')}
            results.append(rec)
    return results


def import_pr115_performance(path: Path) -> list[dict]:
    """Accept a real exported upstream snapshot, keeping raw data and warnings.

    Does not fetch papers, assume clinical independence, or grant semantic trust.
    A digest of the immutable input is stored; unsupported rows are not discarded.
    """
    doc=json.loads(path.read_text())
    raw_rows=doc.get('rows',doc.get('prediction_performance',{}).get('rows'))
    if not isinstance(raw_rows,list):raise ValueError('Expected prediction_performance rows list')
    out=[]
    for index,raw in enumerate(raw_rows):
        if not isinstance(raw,dict):raise ValueError('Each upstream row must be a mapping')
        estimate={'metric':raw.get('metric'),'value':raw.get('value'),'interval':None}
        if raw.get('interval_low') is not None or raw.get('interval_high') is not None:
            estimate['interval']={'lower':raw.get('interval_low'),'upper':raw.get('interval_high'),
                                  'kind':raw.get('interval_type') or 'unknown_reported_interval'}
        context={'endpoint':raw.get('outcome_type'), 'population':raw.get('population'),
                 'country':raw.get('country'),'sex':raw.get('sex'),'age_range':raw.get('age_range'),
                 'horizon_years':None,'horizon_basis':None,'cohort_id':raw.get('cohort'),
                 'validation_set':raw.get('validation_set'),'horizon_original':raw.get('horizon')}
        # Textual horizon requires explicit unit-aware review; no guess from digits.
        out.append({'schema_version':SCHEMA_VERSION,'id':'pr115:'+str(raw.get('id',index)),
          'record_type':'prediction_performance','context':context,
          'measurement':{'id':raw.get('added_component_marker_id'),'assay':None,'matrix':None,'units':None},
          'model':{'id':raw.get('model'),'base_id':raw.get('base_model'),
                   'component_kind':raw.get('added_component_kind'),
                   'added_components':raw.get('component_markers') or [raw.get('added_component')],
                   'upstream_semantic_status':raw.get('checks',{}).get('semantic')},
          'estimate':estimate,'source':{'kind':'literature_extraction','repository':'aleks-kudriavtsev/aging_biomarkers',
             'commit':UPSTREAM_COMMIT,'path':path.name,'locator':f'/rows/{index}',
             'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'pmid':raw.get('source',{}).get('pmid'),
             'review':'upstream_model_semantic_verdict_not_independent_expert_review',
             'synthetic':raw.get('synthetic_source_text',False) or doc.get('not_evidence',False),
             'original':deepcopy(raw.get('source',{}))},'raw_record':deepcopy(raw),
          'clinical_use_ready':False,'usable_as_mortality_coefficient':False})
    return out


def build(root: Path, out: Path, upstream: Path | None=None) -> dict:
    records=literature_records(root)+empirical_records(root)
    if upstream:records.extend(import_pr115_performance(upstream))
    ids=[r['id'] for r in records]
    if len(set(ids))!=len(ids):raise ValueError('Duplicate record identity')
    for row in records:row['audit']=audit_record(row)
    summary={'schema_version':SCHEMA_VERSION,'status':'migration_not_new_clinical_validation',
      'source_commit':BASE_COMMIT,'reviewed_PR115_commit':UPSTREAM_COMMIT,
      'upstream_snapshot_imported':upstream is not None,
      'records':len(records),'literature_records':sum(r['id'].startswith('literature:') for r in records),
      'empirical_records':sum(r['id'].startswith('empirical:') for r in records),
      'retained':sum(r['audit']['retain'] for r in records),
      'records_with_blockers':sum(bool(r['audit']['errors']) for r in records),
      'context_or_source_gaps':sum(bool(r['audit']['gaps']) for r in records),
      'clinical_ready':0,'new_model_fit':False,'new_independent_validation':False}
    out.mkdir(parents=True,exist_ok=True)
    for filename,value in [('evidence_exchange.json',{'summary':summary,'records':records}),
                           ('estimate_contexts.json',project_estimates(records)),('exchange_summary.json',summary)]:
        (out/filename).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    with (out/'evidence_inventory.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['id','kind','source_kind','metric','value','endpoint','sex','age_range','country','horizon_years','interval_kind','errors','gaps'])
        for r in records:
            e,c,a=r['estimate'],r['context'],r['audit']
            # Numeric cells remain numeric, text values cannot become spreadsheet formulas.
            cells=[r['id'],r['record_type'],r['source']['kind'],e['metric'],e['value'],c['endpoint'],c['sex'],str(c['age_range']),c['country'],c['horizon_years'],(e.get('interval') or {}).get('kind'),';'.join(a['errors']),';'.join(a['gaps'])]
            w.writerow(["'"+v if isinstance(v,str) and v.startswith(('=','+','-','@')) else v for v in cells])
    return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--pr115-snapshot',type=Path)
    args=parser.parse_args();print(json.dumps(build(args.root,args.out,args.pr115_snapshot),indent=2))
