"""Isolated synthetic reproductions of reviewed PR115 code paths.

NOT the complete upstream test suite, NOT a real scientific run. The predicates
below reproduce small read source fragments, not an upstream dependency tree.
Source: aging_biomarkers@3c487a8d03b15a49ce1d067e8eba903dab826b10,
 tools/mortality_graph/panel_performance.py (validate_performance),
 tools/mortality_graph/panel_conference_report.py (_added_value, final_panel).
Numeric-gate success is NOT semantic acceptance: upstream semantic review is pending.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import re
from .evidence_contract import METRIC_KINDS, compare_contexts, number_matches


def legacy_numeric_gate(value: float, quote: str) -> bool:
    tokens={float(t) for t in re.findall(r'(?<![\w.])\d+(?:\.\d+)?',quote)}
    return any(math.isclose(abs(value),t,rel_tol=1e-7,abs_tol=1e-9)
        or math.isclose(abs(value)*100,t,rel_tol=1e-7,abs_tol=1e-9)
        or math.isclose(abs(value)*1000,t,rel_tol=1e-7,abs_tol=1e-9) for t in tokens)


def legacy_attachment(performance: dict, marker_id: str, outcome: str, population: str) -> bool:
    # Exactly the selection dimensions of _added_value; rendering omitted.
    clean=lambda value:' '.join(str(value or '').split())
    return (performance.get('added_component_marker_id')==marker_id
       and performance['evidence_kind']=='incremental_value'
       and performance['checks'].get('semantic')=='supported'
       and clean(performance['outcome']).casefold() in {outcome.casefold()}
       and clean(performance['population']).casefold() in {population.casefold()})


def run_probes() -> dict:
    context={'endpoint':'all_cause_mortality','population':'synthetic_population',
      'country':'USA','sex':'female','age_range':[40,79],'horizon_years':10,
      'horizon_basis':'prediction_horizon','cohort_id':'synthetic_A','validation_set':'temporal_validation'}
    p={'added_component_marker_id':'synthetic_marker','evidence_kind':'incremental_value',
       'checks':{'semantic':'supported'},'outcome':'all cause mortality','population':'synthetic_population'}
    rows=[
      {'id':'negative_net_benefit','synthetic_value':.1-.4*.3/.7,
       'upstream_gate_retains':0<=.1-.4*.3/.7<=1,'correct_behavior':'retain negative value with threshold and action'},
      {'id':'NRI_above_one','synthetic_value':1.2,'upstream_gate_retains':-1<=1.2<=1,
       'correct_behavior':'retain total NRI within [-2,2]'},
      {'id':'sign_flip_numeric_gate','quote':'delta C = -0.01','extracted_value':.01,
       'upstream_numeric_gate_passes':legacy_numeric_gate(.01,'delta C = -0.01'),
       'new_numeric_gate_passes':number_matches(.01,'-0.01')},
      {'id':'undeclared_scale_numeric_gate','quote':'delta C = 4','extracted_value':.004,
       'upstream_numeric_gate_passes':legacy_numeric_gate(.004,'delta C = 4'),
       'new_numeric_gate_passes':number_matches(.004,'4')},
      {'id':'Brier_taxonomy','upstream_kind':'calibration','new_kind':METRIC_KINDS['Brier_score']},
    ]
    for axis,value in [('horizon_years',1),('sex','male'),('cohort_id','synthetic_B')]:
        changed={**context,axis:value}
        rows.append({'id':'attachment_mismatch_'+axis,
          'upstream_attachment_selects':legacy_attachment({**p,**changed,'outcome':p['outcome']},'synthetic_marker',p['outcome'],p['population']),
          'new_context_status':compare_contexts({'context':context},{'context':changed})['status']})
    estimates=[{'id':'synthetic_M','measure':'HR','value':1.2,'sex':'male','age_range':[40,59],'country':'USA'},
               {'id':'synthetic_F','measure':'HR','value':1.6,'sex':'female','age_range':[60,79],'country':'DEU'}]
    keys=('id','measure','value','ci_low','ci_high','contrast','horizon','followup_basis','followup_years','cohort','study_design')
    projected=[{k:c.get(k) for k in keys} for c in estimates]
    rows.append({'id':'per_estimate_applicability_projection','upstream_estimate_projection_preserves_sex':all('sex' in r for r in projected),
      'upstream_raw_evidence_still_has_context':True,'new_requirement':'sex, age and country must accompany each effect'})
    assert len(rows)==9
    assert rows[0]['upstream_gate_retains'] is False and rows[1]['upstream_gate_retains'] is False
    assert rows[2]['upstream_numeric_gate_passes'] and not rows[2]['new_numeric_gate_passes']
    assert rows[3]['upstream_numeric_gate_passes'] and not rows[3]['new_numeric_gate_passes']
    assert all(r['upstream_attachment_selects'] and r['new_context_status']=='incompatible' for r in rows[5:8])
    return {'status':'isolated_predicate_reproductions_not_full_upstream_run','synthetic_only':True,
      'not_clinical_evidence':True,'reviewed_commit':'3c487a8d03b15a49ce1d067e8eba903dab826b10',
      'source_blob_ids':{'panel_performance.py':'4d0ceb4aa4769eecdbc6a3922807aa6ee5114660',
                         'panel_conference_report.py':'67ed637c0fa491e9321fd3e4652afe6944aabc63'},
      'important_limit':'numeric gate success does not establish upstream semantic acceptance',
      'probes':rows}

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(run_probes(),indent=2,ensure_ascii=False)+'\n')
    print('9 isolated synthetic probes completed; not an end-to-end PR115 run')
