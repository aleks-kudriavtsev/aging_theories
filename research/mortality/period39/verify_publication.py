#!/usr/bin/env python3
"""Isolated publication audit: unchanged original arithmetic; explicit published Git provenance.

Historical verify_results.py remains unchanged. No historical ancestor is fabricated.

Reads only already-authorized local public population aggregates, pinned docs,
output CSV/JSON, chart bytes and saved logs. No network or person records.
"""
import argparse
import csv
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import re
import subprocess
import sys

PROTOCOL='c0a570e0ac40de94da8a05b53f5e9d80aaf7aacf'
SCIENCE='151fe91b1540b91926baf3b3cc782cc7515ae55b'
SOURCE='31fff9249237a33d9a10ac107d6710a04bfec517'
PREFIX='research/mortality/period39/'
CODE_FILES={PREFIX+x for x in ('PROTOCOL_RU.md','STUDY_PLAN.json','SOURCE_MANIFEST.json','METHOD_REVIEW.md','SHFM_REVIEW.md','SHFM_RIGHTS.md','REVIEW_MANIFEST.json','IMPLEMENTATION_NOTES_RU.md','__init__.py','reference.py','analyze.py','plot.py')}|{'tests/test_period39.py','tests/test_analyze39.py'}
SOURCE_FILES={'data/mortality/mx_qx_usa_deu_rus.csv','data/mortality/README.md','docs/12_cross_country_mortality.md','tools/fetch_wpp_mortality.py'}
CSV_FIELDS=['country','source_year','source_sex','age_exact_years','horizon_years','end_exact_birthday','CF_probability','CF_rounding_lower','CF_rounding_upper','UDD_ax_0_5_probability','UDD_rounding_lower','UDD_rounding_upper','UDD_minus_CF','stored_author_qx_product','stored_qx_minus_recomputed_CF','output_scope','clinical_validation','current_calendar_forecast']
METHODS=('constant_force','UDD_ax_0_5')
HALF=Decimal('0.0000005')
DECIMAL_TOL=Decimal('3e-83')


def sha(raw):return hashlib.sha256(raw).hexdigest()
def Q(x):return x if type(x)is Fraction else Fraction(str(x))
def loads(raw):
    def pairs(items):
        d={}
        for k,v in items:
            if k in d:raise ValueError('Duplicate JSON key '+k)
            d[k]=v
        return d
    return json.loads(raw,object_pairs_hook=pairs,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite JSON '+x)))
def blob(repo,commit,name):return subprocess.check_output(['git','-C',str(repo),'show',commit+':'+name],stderr=subprocess.PIPE)

class Audit:
    def __init__(self):self.total=0;self.failures=[];self.max_decimal_error=Decimal(0);self.max_float_error=0.
    def check(self,name,ok):
        self.total+=1
        if not ok:self.failures.append({'id':name})
    def equal(self,name,a,b):self.check(name,type(a)==type(b) and a==b)
    def close(self,name,x,y,tol=2e-15):
        valid=type(x)in(int,float) and math.isfinite(x);error=abs(x-y) if valid else math.inf
        if math.isfinite(error):self.max_float_error=max(self.max_float_error,error)
        self.check(name,valid and error<=tol)
    def decimal(self,name,observed,expected):
        try:
            val=Decimal(observed)
            with localcontext() as c:c.prec=95;error=abs(val-expected)
            valid=val.is_finite() and error<=DECIMAL_TOL
            if error.is_finite():self.max_decimal_error=max(self.max_decimal_error,error)
        except Exception:valid=False
        self.check(name,valid)


def exact_product(values):
    result=Fraction(1)
    for value in values:result*=value
    return result

def oracle(mx,qx):
    """Independent CF via annual exp products; UDD/qx via exact rational products."""
    with localcontext() as c:
        c.prec=85
        lower=[max(Decimal(0),m-HALF) for m in mx];upper=[m+HALF for m in mx]
        cf=[];udd=[]
        for rates in (mx,lower,upper):
            surv=Decimal(1)
            for rate in rates:surv*=(-rate).exp()
            cf.append(1-surv)
            us=exact_product((2-Q(m))/(2+Q(m)) for m in rates)
            udd.append(Decimal(1)-Decimal(us.numerator)/Decimal(us.denominator))
        qs=exact_product(1-Q(p) for p in qx)
        diagnostic=Decimal(1)-Decimal(qs.numerator)/Decimal(qs.denominator)
        return {'constant_force':dict(zip(('probability','lower','upper'),cf)), 'UDD_ax_0_5':dict(zip(('probability','lower','upper'),udd)), 'stored_qx_product':+diagnostic,'delta_vs_constant_force':+(diagnostic-cf[0])}


def audit_results(repo,results,check_charts=True,code_ref="HEAD"):
    a=Audit();root=repo/PREFIX
    published=subprocess.check_output(['git','-C',str(repo),'rev-parse',code_ref],text=True).strip()
    cmraw=(root/'CODE_MANIFEST.json').read_bytes();cm=loads(cmraw)
    a.equal('publication.original_manifest_SHA',sha(cmraw),'f1aa57c37244fdbf06def5e2361696fabc7b4bdd63ed4159be9eeaada64bf8f7')
    a.equal('code.schema',cm['schema'],'period39-code-manifest/1');a.equal('code.protocol',cm['protocol_commit'],PROTOCOL)
    a.equal('code.inventory14',set(x['path'] for x in cm['files']),CODE_FILES);a.equal('code.count14',len(cm['files']),14)
    for item in cm['files']:
        name=item['path'];a.equal('code.current.'+name,sha((repo/name).read_bytes()),item['sha256']);a.equal('code.published_Git.'+name,sha(blob(repo,published,name)),item['sha256'])
    protocol_names=[PREFIX+x for x in ('PROTOCOL_RU.md','STUDY_PLAN.json','SOURCE_MANIFEST.json')]
    for name in protocol_names:a.equal('protocol.published_unchanged.'+name,sha(blob(repo,published,name)),next(x['sha256'] for x in cm['files'] if x['path']==name))
    planraw=(root/'STUDY_PLAN.json').read_bytes();plan=loads(planraw);smraw=(root/'SOURCE_MANIFEST.json').read_bytes();sm=loads(smraw)
    a.equal('source.manifest_plan_pin',sha(smraw),plan['source']['manifest_sha256']);a.equal('source.commit',sm['source_commit'],SOURCE);a.equal('source.exact_files4',set(x['path'] for x in sm['files']),SOURCE_FILES);a.equal('source.count4',len(sm['files']),4)
    cached={};source_root=root/'source_snapshot'
    for item in sm['files']:
        name=item['path'];cached[name]=(source_root/name).read_bytes();a.equal('source.current.'+name,sha(cached[name]),item['sha256']);a.equal('source.published_snapshot_Git.'+name,sha(blob(repo,published,PREFIX+'source_snapshot/'+name)),item['sha256'])
    a.equal('source.original_RDA_hash_unknown',sm['official_raw_RDA_sha256'],None);a.equal('source.no_new_downloads',sm['new_raw_source_downloads'],False);a.equal('source.qx_author_derived',sm['qx_is_author_derived_not_official'],True);a.equal('source.declared_unit',sm['declared_mx_unit'],'annual_deaths_per_person_year');a.equal('source.unit_evidence',sm['verified_unit_level'],'local_extractor_and_repo_docs; official_raw_RDA_numeric_scaling_not_reverified')
    reviews=loads((root/'REVIEW_MANIFEST.json').read_bytes());a.equal('reviews.schema',reviews['schema'],'period39-review-manifest/1');a.equal('reviews.files3',{x['path'] for x in reviews['files']},{'METHOD_REVIEW.md','SHFM_REVIEW.md','SHFM_RIGHTS.md'})
    for x in reviews['files']:a.equal('reviews.sha.'+x['path'],sha((root/x['path']).read_bytes()),x['sha256'])
    a.equal('reviews.no_RDA_bytes',reviews['original_RDA_bytes_retrieved'],False);a.equal('reviews.no_SHFM5_reproduction',reviews['SHFM_5year_model_reproduced'],False)
    shfm=(root/'SHFM_REVIEW.md').read_text();rights=(root/'SHFM_RIGHTS.md').read_text();a.check('reviews.SHFM_updated_2013_unqualified','2013' in shfm and 'not yet qualified' in shfm and 'Updated2–5year baselines | Blocked' in shfm);a.check('reviews.SHFM_licence_not_granted','BLOCKED' in rights and 'no general open-source' in rights.lower())
    if a.failures:raise ValueError('Pre-parse source/freeze failure: '+str([x['id'] for x in a.failures]))
    reader=csv.DictReader(io.StringIO(cached[plan['source']['csv']].decode()));a.equal('source.CSV_schema',reader.fieldnames,['country','sex','year','age','mx','qx','series']);target={}
    for row in reader:
        if row['country']=='RUS' and row['year']=='2023' and 40<=int(row['age'])<=89:
            key=(row['sex'],int(row['age']))
            if key in target:raise ValueError('Duplicate target source key')
            if row['sex'] not in ('F','M') or row['series']!='estimate':raise ValueError('Changed target source category')
            mx,qx=Decimal(row['mx']),Decimal(row['qx']);target[key]=(mx,qx)
    a.equal('source.target100',len(target),100);a.equal('source.exact_target_keys',set(target),set(itertools.product(('F','M'),range(40,90))))
    max_qx_year_error=Decimal(0)
    for key,(mx,qx) in target.items():
        tag='source.target.'+str(key);a.check(tag+'.positive_valid_UDD',mx>0 and mx+HALF<=2);a.check(tag+'.6decimal',mx.as_tuple().exponent>=-6 and qx.as_tuple().exponent>=-6);a.check(tag+'.qx_unit_interval',0<=qx<=1)
        with localcontext() as c:c.prec=85;error=abs(qx-(1-(-mx).exp()))
        max_qx_year_error=max(max_qx_year_error,error);a.check(tag+'.author_qx_rounding_diagnostic',error<=Decimal('0.000001'))
    nested=results/'aggregate/RESULTS_SUMMARY.json';flat=results/'RESULTS_SUMMARY.json'
    if nested.exists() and flat.exists() and nested.read_bytes()!=flat.read_bytes():raise ValueError('Ambiguous nested and flat JSON bytes differ')
    output=nested if nested.exists() else flat;csvfile=output.parent/'REFERENCE_TABLE.csv';raw=output.read_bytes();csvraw=csvfile.read_bytes();r=loads(raw)
    a.equal('results.schema',r['schema'],'period39-results/1');a.equal('results.protocol',r['protocol_commit'],PROTOCOL);a.equal('results.plan_sha',r['plan_sha256'],sha(planraw));a.equal('results.code_manifest_sha',r['code_manifest_sha256'],sha(cmraw));a.equal('results.CSV_sha',r['reference_csv_sha256'],sha(csvraw))
    a.equal('target.country',plan['target']['country'],'RUS');a.equal('target.ages40_80',plan['target']['exact_baseline_ages_inclusive'],[40,80]);a.equal('target.sexes',plan['target']['source_sexes'],['F','M']);a.equal('target.horizons',plan['target']['horizons_years'],[5,10]);a.equal('results.source_rows100',r['source_rows'],100);a.equal('results.period_year',r['source_year'],2023);a.equal('results.period_assumption',r['period_assumption'],'hold_2023_attained_age_schedule');a.equal('results.scope',r['scope'],'provisional_local_derived_period_reference_under_declared_units');a.equal('results.mixed_population',r['source_population'],'allcause_mixed_general_population_including_disease')
    for flag in ('rounding_bounds_are_confidence_intervals','model_assumption_difference_is_uncertainty_CI','disease_free_baseline','original_RDA_numeric_scale_reverified','current_calendar_forecast','biomarker_adjusted','clinical_validation','models_fitted','SHFM_updated_5year_reproduced'):a.equal('results.false.'+flag,r[flag],False)
    a.equal('results.no_new_person_inputs',r['new_person_level_inputs'],0)
    expected_keys=list(itertools.product(('F','M'),range(40,81),(5,10)));rows=r['rows'];a.equal('results.rows164',len(rows),164);a.equal('results.row_count_claim',r['reference_rows'],164)
    hashes={x['path']:x['sha256'] for x in sm['files']};bykey={};oracles={};max_width={m:0. for m in METHODS};max_difference=0.;max_stored_delta=0.
    expected_metadata={'research_only':True,'reference_status':'provisional_local_derived_population_period_reference','clinical_forecast':False,'clinical_validation':False,'individual_calibration_performed':False,'source_population_includes_diseased':True,'disease_free_or_stage_specific_baseline':False,'forecast_of_future_calendar_rates':False,'future_schedule_held_at_source2023':True,'source_rates_declared_unit':'annual_deaths_per_person_year','source_unit_evidence':sm['verified_unit_level'],'official_raw_RDA_numeric_scaling_reverified':False,'official_raw_RDA_checksum_verified':False,'source_files_verified':True,'source_commit':SOURCE,'source_hashes':hashes,'stored_qx_is_author_derived_not_official_UN_qx':True,'rounding_bounds_are_confidence_intervals':False,'rounding_assumption':'each_mx_rounded_to_nearest6decimals_with_independent_closed_halfwidth_0.5e-6','assumption_difference_is_model_selection':False,'interpolation_or_extrapolation_performed':False,'sexes_averaged':False,'new_outcomes_or_participant_records_read':False}
    for i,row in enumerate(rows):
        req=row['request'];key=(req['source_sex'],req['age_exact_years'],req['horizon_years']);tag='row.'+str(i);a.equal(tag+'.order_and_key',key,expected_keys[i]);a.check(tag+'.no_duplicate',key not in bykey);bykey[key]=row;s,age,h=key
        a.equal(tag+'.exact_row_schema',set(row),{'methods','diagnostic','precision','years_count','schema','request','source_ages','end_exact_birthday','metadata'})
        request={'schema':'period39-request/1','model':'RUS_WPP2023_period_reference','country':'RUS','source_year':2023,'source_sex':s,'age_exact_years':age,'horizon_years':h,'age_origin':'exact_birthday','period_assumption':'hold_2023_attained_age_schedule','purpose':'population_period_reference_only'}
        a.equal(tag+'.exact_request',req,request);a.check(tag+'.integer_coordinates',type(age)is int and type(h)is int and type(req['source_year'])is int);ages=list(range(age,age+h));a.equal(tag+'.attained_ages',row['source_ages'],ages);a.equal(tag+'.end_birthday',row['end_exact_birthday'],age+h);a.equal(tag+'.year_count',row['years_count'],h);a.equal(tag+'.metadata',row['metadata'],expected_metadata);a.equal(tag+'.precision',row['precision'],{'decimal_digits':85,'mx_rounding_halfwidth':'0.0000005'})
        values=[target[s,x] for x in ages];truth=oracle([x[0] for x in values],[x[1] for x in values]);oracles[key]=truth;a.equal(tag+'.method_names',set(row['methods']),set(METHODS))
        for method in METHODS:
            o=row['methods'][method];ref=truth[method];a.equal(tag+'.method_schema.'+method,set(o),{'probability','rounding_bounds','decimal'});a.equal(tag+'.bounds_schema.'+method,set(o['rounding_bounds']),{'lower','upper'})
            for field in ('probability','lower','upper'):
                a.decimal(tag+'.85digit.'+method+'.'+field,o['decimal'][field],ref[field]);fl=o['probability'] if field=='probability' else o['rounding_bounds'][field];a.close(tag+'.float.'+method+'.'+field,fl,float(ref[field]))
            lo,p,hi=(Decimal(o['decimal'][x]) for x in ('lower','probability','upper'));a.check(tag+'.bounds_order.'+method,0<=lo<=p<=hi<=1);a.check(tag+'.width_guard.'+method,hi-lo<=Decimal(h)*Decimal('0.000001')*Decimal('1.00001'));width=o['rounding_bounds']['upper']-o['rounding_bounds']['lower'];max_width[method]=max(max_width[method],width)
        a.check(tag+'.UDD_not_less_than_CF',truth['UDD_ax_0_5']['probability']>=truth['constant_force']['probability'])
        for field in ('stored_qx_product','delta_vs_constant_force'):a.decimal(tag+'.diagnostic_decimal.'+field,row['diagnostic']['decimal'][field],truth[field]);a.close(tag+'.diagnostic_float.'+field,row['diagnostic'][field],float(truth[field]))
        max_difference=max(max_difference,abs(row['methods']['UDD_ax_0_5']['probability']-row['methods']['constant_force']['probability']));max_stored_delta=max(max_stored_delta,abs(row['diagnostic']['delta_vs_constant_force']))
    a.equal('results.exact164keys',set(bykey),set(expected_keys));a.equal('results.selected_rows20',len(r['selected_rows']),20);a.equal('results.selected_rows_exact',r['selected_rows'],[row for row in rows if row['request']['age_exact_years'] in (40,50,60,70,80)])
    for sex in ('F','M'):
        a.equal('age80_10_support.'+sex,bykey[sex,80,10]['source_ages'],list(range(80,90)));a.equal('age80_10_endpoint.'+sex,bykey[sex,80,10]['end_exact_birthday'],90)
    chains=r['chaining_rows'];a.equal('chains.rows72',len(chains),72);a.equal('chains.count_claim',r['chaining_age_sex_rows'],72);a.equal('chains.method_count144',r['chaining_scalar_method_checks'],144);a.equal('horizon.count164',r['horizon_coherence_scalar_checks'],164)
    chain_keys=list(itertools.product(('F','M'),range(40,76)))
    max_chain_error=0.
    with localcontext() as c:
        c.prec=85
        for i,row in enumerate(chains):
            s,age=row['source_sex'],row['age_exact_years'];a.equal('chain.'+str(i)+'.key',(s,age),chain_keys[i]);a.equal('chain.'+str(i)+'.methods',set(row['methods']),set(METHODS))
            for method in METHODS:
                pa,pb,total=(oracles[k][method]['probability'] for k in ((s,age,5),(s,age+5,5),(s,age,10)));a.check('chain.'+str(i)+'.85digit.'+method,abs(total-(pa+(1-pa)*pb))<=Decimal('1e-82'))
                x,y,z=(bykey[k]['methods'][method]['probability'] for k in ((s,age,5),(s,age+5,5),(s,age,10)));error=abs(z-(x+(1-x)*y));max_chain_error=max(max_chain_error,error);a.close('chain.'+str(i)+'.float_error.'+method,row['methods'][method]['absolute_error'],error,tol=1e-30);a.equal('chain.'+str(i)+'.declared.'+method,row['methods'][method]['passed'],True)
        for s,age,method in itertools.product(('F','M'),range(40,81),METHODS):a.check('horizon.'+s+str(age)+method,oracles[s,age,10][method]['probability']>=oracles[s,age,5][method]['probability'])
    declines=[]
    for s,h,age in itertools.product(('F','M'),(5,10),range(40,80)):
        if bykey[s,age+1,h]['methods']['constant_force']['probability']<bykey[s,age,h]['methods']['constant_force']['probability']:declines.append({'source_sex':s,'age_from':age,'horizon_years':h})
    a.equal('diagnostics.age_decreases',r['age_decreases_diagnostic_only'],declines)
    for method in METHODS:a.close('diagnostics.max_width.'+method,r['rounding_envelope_max_width'][method],max_width[method],tol=1e-30)
    a.close('diagnostics.max_method_difference',r['max_UDD_CF_absolute_difference'],max_difference,tol=1e-30);a.close('diagnostics.max_stored_qx_difference',r['max_stored_qx_CF_absolute_difference'],max_stored_delta,tol=1e-30)
    table_reader=csv.DictReader(io.StringIO(csvraw.decode()));a.equal('output_CSV.exact_schema',table_reader.fieldnames,CSV_FIELDS);table=list(table_reader);a.equal('output_CSV.rows164',len(table),164)
    for i,(csvrow,row) in enumerate(zip(table,rows)):
        s,age,h=expected_keys[i];tag='CSV.'+str(i);a.equal(tag+'.key',(csvrow['source_sex'],int(csvrow['age_exact_years']),int(csvrow['horizon_years'])),(s,age,h));a.equal(tag+'.context',[csvrow['country'],csvrow['source_year'],csvrow['end_exact_birthday'],csvrow['output_scope'],csvrow['clinical_validation'],csvrow['current_calendar_forecast']],['RUS','2023',str(age+h),'provisional_local_derived_period_reference','False','False'])
        for method,prefix in [('constant_force','CF'),('UDD_ax_0_5','UDD')]:
            names={'probability':'CF_probability' if method=='constant_force' else 'UDD_ax_0_5_probability','lower':prefix+'_rounding_lower','upper':prefix+'_rounding_upper'}
            for field,col in names.items():a.equal(tag+'.decimal_exact_JSON.'+col,csvrow[col],row['methods'][method]['decimal'][field]);a.decimal(tag+'.independent.'+col,csvrow[col],oracles[s,age,h][method][field])
        for col,field in [('stored_author_qx_product','stored_qx_product'),('stored_qx_minus_recomputed_CF','delta_vs_constant_force')]:a.equal(tag+'.diagnostic_exact.'+col,csvrow[col],row['diagnostic']['decimal'][field]);a.decimal(tag+'.diagnostic_oracle.'+col,csvrow[col],oracles[s,age,h][field])
        a.close(tag+'.method_difference',float(csvrow['UDD_minus_CF']),row['methods']['UDD_ax_0_5']['probability']-row['methods']['constant_force']['probability'],tol=1e-30)
    log=loads((results/'run.log').read_bytes());a.equal('log.completed',log['status'],'completed_local_period_reference');a.equal('log.output_SHA',log['result_sha256'],sha(raw));a.equal('log.CSV_SHA',log['csv_sha256'],sha(csvraw));a.equal('log.rows164',log['reference_rows'],164);a.equal('log.not_clinical',log['clinical_validation'],False)
    for v in ('312','313'):
        text=(results/('tests_python'+v+'.log')).read_text();a.check('testlog.label'+v+'.27_completed',bool(re.search(r'Ran 27 tests\b',text)));a.check('testlog.label'+v+'.OK',bool(re.search(r'^OK\s*$',text,re.M)))
    chart_hashes={}
    if check_charts:
        logpath=results/'chart_clean.log'
        if not logpath.exists():logpath=results/'chart.log'
        candidates=[]
        for line in logpath.read_text().splitlines():
            try:doc=loads(line)
            except Exception:continue
            if isinstance(doc,dict) and doc.get('schema')=='period39-chart-provenance/1':candidates.append(doc)
        if len(candidates)!=1:raise ValueError('Exactly one chart provenance JSON record required')
        chart=candidates[0];a.equal('chart.source_result_SHA',chart['source_results_sha256'],sha(raw));a.equal('chart.not_clinical',chart['clinical_validation'],False)
        copied_provenance=results/'CHART_PROVENANCE.json'
        a.check('chart.optional_saved_provenance_matches_clean_log',not copied_provenance.exists() or all(loads(copied_provenance.read_bytes()).get(k)==v for k,v in chart.items()))
        for ext,key in [('svg','svg_sha256'),('png','png_sha256')]:
            paths=[results/('PERIOD_REFERENCE.'+ext),results/('PERIOD_REFERENCE_CACHECHECK.'+ext)];present=[p for p in paths if p.exists()]
            if not present:raise ValueError('Missing chart '+ext)
            a.check('chart.all_present_file_SHAs.'+ext,all(sha(p.read_bytes())==chart[key] for p in present))
            chart_hashes[key]=chart[key];data=present[0].read_bytes();a.check('chart.inline_result_SHA.'+ext,sha(raw).encode() in data)
            if ext=='svg':a.check('chart.SVG_signature',b'<svg' in data and b'not individual prognosis or clinical validation' in data)
            else:a.check('chart.PNG_signature',data.startswith(b'\x89PNG\r\n\x1a\n'))
    forbidden={'participant_id','person_id','SEQN','biomarker_values','participant_records','disease_specific_risk','individual_clinical_risk'}
    def find_forbidden(x):
        if isinstance(x,dict):return [key for key in x if key in forbidden]+sum((find_forbidden(v) for v in x.values()),[])
        if isinstance(x,list):return sum((find_forbidden(v) for v in x),[])
        return []
    a.equal('export.no_personal_or_disease_risk_fields',find_forbidden(r),[])
    return {'schema':'period39-independent-compact-audit/1','passed':not a.failures,'checks_passed':a.total-len(a.failures),'checks_total':a.total,'failures':a.failures,'scientific_commit':SCIENCE,'protocol_commit':PROTOCOL,'source_commit':SOURCE,'publication_code_ref':published,'historical_Git_ancestry_revalidated':False,'execution_origin':'local','producer_imports_or_replays':False,'new_network_or_person_data':False,'charts_checked':check_charts,'hashes':{'result_sha256':sha(raw),'reference_csv_sha256':sha(csvraw),'code_manifest_sha256':sha(cmraw),'source_manifest_sha256':sha(smraw),'plan_sha256':sha(planraw),**chart_hashes},'checked_counts':{'code_pins':14,'protocol_files':3,'source_pins':4,'target_source_rows':100,'reference_cells':164,'selected_reference_cells':20,'methods_per_cell':2,'chaining_rows':72,'chaining_method_checks':144,'horizon_coherence_checks':164,'review_pins':3},'independent_reference':{'decimal_digits':85,'CF_form':'product_of_annual_exponential_survivals','UDD_form':'exact_Fraction_product_(2-m)/(2+m)','stored_qx_form':'exact_Fraction_survival_product','max_absolute_decimal_error':str(a.max_decimal_error),'decimal_absolute_comparison_tolerance':str(DECIMAL_TOL),'max_float_error':a.max_float_error,'max_one_year_stored_qx_rounding_error':str(max_qx_year_error),'max_float_chain_error':max_chain_error},'diagnostics':{'max_UDD_CF_absolute_difference':max_difference,'max_stored_qx_CF_absolute_difference':max_stored_delta,'rounding_envelope_max_width':max_width,'age_decreases_count':len(declines)},'limits':['Conditional provisional local-derived population period reference under declared annual units; original UN RDA bytes and numeric scaling are not reauthenticated. Hashes authenticate published cached files and original manifest; historical Git ancestry is not revalidated in this isolated publication.','Calendar2023 rates are held at future attained ages. This is not a current-calendar or cohort forecast, disease-free/stage baseline, personalized prediction, biomarker model or clinical validation. Mixed-population all-cause mortality must not be added to postdiagnosis all-cause risk.','CF and UDD are predeclared within-age assumptions. Rounded-mx box envelopes are arithmetic sensitivity bounds, not stochastic confidence intervals or forecasts; method differences do not quantify all model uncertainty.','SHFM updated2013 five-year reproduction remains blocked/unperformed; source/rights reviews are preserved, with no model code or numerical baselines inferred.','Saved local logs labelled Python312/313 show27 tests/OK; plain log text alone does not authenticate interpreter binaries or remote CI. Chart checks authenticate bytes/source-SHA linkage, not independent rendering or epidemiological validity.']}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,required=True);p.add_argument('--results',type=Path,required=True,help='Run folder with aggregate/RESULTS_SUMMARY.json + CSV, or flat saved-results folder; run/test/chart logs and plot files alongside');p.add_argument('--out',type=Path,required=True);p.add_argument('--code-ref',default='HEAD');x=p.parse_args();r=audit_results(x.repo.resolve(),x.results.resolve(),code_ref=x.code_ref);x.out.parent.mkdir(parents=True,exist_ok=True);x.out.write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    lines=['# Period39 independent compact audit','',f"**{'PASS' if r['passed'] else 'FAIL'}: {r['checks_passed']}/{r['checks_total']} checks.**",'',f"Scientific code `{SCIENCE}`; protocol `{PROTOCOL}`.",'','Checked14 current/Git code pins,3 hash-preserved protocol files; historical ancestry is not revalidated,4 cached aggregate/source-context pins against published Git,100 target source keys, all164 reference cells/CSV rows,72 chaining rows×2methods and164 horizon checks. Three source-review hashes and blocked SHFM status are preserved.','',f"Independent85-digit Decimal/Fraction oracle maximum decimal error: `{r['independent_reference']['max_absolute_decimal_error']}`.",'','No producer imports/replays, network, source retries, fitting or participant processing.','', 'Limits:','']+['- '+s for s in r['limits']]+['','Hashes:','']+['- `'+k+'`: `'+v+'`' for k,v in r['hashes'].items()]
    if r['failures']:lines+=['','Failed checks:','']+['- `'+v['id'].replace('|','\\|')+'`' for v in r['failures']]
    x.out.with_suffix('.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'passed':r['passed'],'checks_passed':r['checks_passed'],'checks_total':r['checks_total'],'out':str(x.out)}));return 0 if r['passed'] else 1

if __name__=='__main__':sys.exit(main())
