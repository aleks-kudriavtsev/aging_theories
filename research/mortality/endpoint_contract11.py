"""Observable endpoint contract, not a fitted mortality model.

Public NCHS LMF2019 dictionary (28 April 2022), pp1-2; description p3.
Codes refer to INITIAL UNDERLYING cause, never contributing-cause flags.
A known-dead person with missing cause remains a death. An unavailable subtype
is not a negative label. No attempt to reconstruct perturbed/restricted records.
"""
from __future__ import annotations
from collections.abc import Sequence
import math
import re

PUBLIC_CAUSES = {
 '001': {'id':'heart','label_ru':'Болезни сердца','ucod113':'054-068',
         'icd10':'I00-I09,I11,I13,I20-I51'},
 '002': {'id':'malignant','label_ru':'Злокачественные новообразования','ucod113':'019-043','icd10':'C00-C97'},
 '003': {'id':'clrd','label_ru':'Хронические болезни нижних дыхательных путей','ucod113':'082-086','icd10':'J40-J47'},
 '004': {'id':'accidents','label_ru':'Несчастные случаи','ucod113':'112-123','icd10':'V01-X59,Y85-Y86'},
 '005': {'id':'cerebrovascular','label_ru':'Цереброваскулярные болезни','ucod113':'070','icd10':'I60-I69'},
 '006': {'id':'alzheimer','label_ru':'Болезнь Альцгеймера','ucod113':'052','icd10':'G30'},
 '007': {'id':'diabetes','label_ru':'Сахарный диабет','ucod113':'046','icd10':'E10-E14'},
 '008': {'id':'influenza_pneumonia','label_ru':'Грипп и пневмония','ucod113':'076-078','icd10':'J10-J18'},
 '009': {'id':'renal_broad','label_ru':'Нефрит, нефротический синдром и нефроз','ucod113':'097-101','icd10':'N00-N07,N17-N19,N25-N27'},
 '010': {'id':'other','label_ru':'Все остальные причины','ucod113':'residual','icd10':'residual'},
 'UNK': {'id':'unknown','label_ru':'Смерть с неизвестной причиной','ucod113':None,'icd10':None},
}
# These are definitions of desired outcomes, NOT trained models or diagnoses.
TARGETS = {
 'ihd': ('001','proper_subset','I20-I25'),
 'heart': ('001','same_group','I00-I09,I11,I13,I20-I51'),
 'malignant': ('002','same_group','C00-C97'),
 'lung_cancer': ('002','proper_subset','C33-C34'),
 'colorectal_cancer': ('002','proper_subset','C18-C21'),
 'clrd': ('003','same_group','J40-J47'),
 'copd_J44': ('003','proper_subset','J44'),
 'accidents': ('004','same_group','V01-X59,Y85-Y86'),
 'all_external': (None,'not_identifiable','V01-Y89'),
 'cerebrovascular': ('005','same_group','I60-I69'),
 'ischemic_stroke': ('005','proper_subset','I63'),
 'alzheimer': ('006','same_group','G30'),
 'dementia_F01_F03': ('010','proper_subset','F01-F03'),
 'diabetes': ('007','same_group','E10-E14'),
 'influenza_pneumonia': ('008','same_group','J10-J18'),
 'renal_broad': ('009','same_group','N00-N07,N17-N19,N25-N27'),
 'ckd_N18': ('009','proper_subset','N18'),
 'noninfectious_natural': (None,'not_identifiable',None),
}


def allowed_codes(survey_start: int) -> frozenset[str]:
    if type(survey_start) is not int or survey_start not in range(1999,2018,2):
        raise ValueError('Specify a supported NHANES start year, 1999..2017 odd years')
    return frozenset({'001','002','010'} if survey_start>=2015 else PUBLIC_CAUSES.keys()-{'UNK'})


def normalize_public_cause(value) -> str | None:
    if value is None or (type(value) is float and math.isnan(value)):
        return None
    if isinstance(value,str):
        text=value.strip()
        if text in {'','.','UNK'}: return None
        if not re.fullmatch(r'\d{3}',text): raise ValueError('Expected three-digit UCOD_LEADING')
        code=text
    elif type(value) in (int,float) and math.isfinite(value) and value==int(value):
        code=f'{int(value):03d}'
    else: raise ValueError('Invalid public cause type')
    if code not in PUBLIC_CAUSES or code=='UNK':raise ValueError('Undeclared public cause code')
    return code


def classify_public(status, cause, *, survey_start: int) -> str:
    """Return alive or one of ten underlying causes or UNK; never censor a death."""
    permitted=allowed_codes(survey_start)
    if type(status) not in (int,float) or status not in (0,1):
        raise ValueError('Known binary mortality status required; eligibility is separate')
    code=normalize_public_cause(cause)
    if status==0:
        if code is not None:raise ValueError('Living status conflicts with cause of death')
        return 'alive'
    if code is None:return 'UNK'
    if code not in permitted:raise ValueError('Cause is unavailable in this survey release')
    return code


def observability(target: str, *, survey_start: int) -> dict:
    if target not in TARGETS:raise ValueError('Unknown target definition')
    code,relation,icd=TARGETS[target]
    observable=relation=='same_group' and code in allowed_codes(survey_start)
    reason=('observable_public_group' if observable else 'source_coarsened_for_cycle'
            if code is not None and code not in allowed_codes(survey_start)
            else 'subtype_not_identified' if relation=='proper_subset' else 'mixed_or_multiple_groups')
    return {'target':target,'icd10_definition':icd,'public_parent':code,'relation':relation,
            'label_observable':observable,'reason':reason,'clinical_use_ready':False,
            'cause_perturbation_possible':True}


def fixed_horizon_states(statuses: Sequence,causes: Sequence,months: Sequence,
                         *,survey_start: int,horizon_years: int) -> list[str]:
    """Empirical complete-ascertainment labels; survivors after a horizon are alive.

    Refuse early censoring instead of coding a short-followed survivor negative.
    No finite-risk extrapolation. This is not an individual prediction routine.
    """
    if type(horizon_years) is not int or horizon_years<=0:raise ValueError('Positive integer horizon required')
    if not (len(statuses)==len(causes)==len(months)) or len(statuses)==0:
        raise ValueError('Nonempty aligned arrays required')
    result=[]
    for status,cause,t in zip(statuses,causes,months):
        state=classify_public(status,cause,survey_start=survey_start)
        if type(t) not in (int,float) or not math.isfinite(t) or t<0:
            raise ValueError('Nonnegative finite follow-up months required')
        if status==0 and t<12*horizon_years:raise ValueError('Early censoring: use survival estimator')
        result.append(state if status==1 and t<=12*horizon_years else 'alive')
    return result


def check_linkage_lineage(mortality_releases: Sequence[str]) -> dict:
    """Scientific input check reflecting CDC RDC public/restricted-LMF prohibition.

    Not an access-control or authorisation mechanism; no restricted files opened.
    Ordinary public NHANES baseline measurements are not public mortality linkage.
    """
    if isinstance(mortality_releases,(str,bytes)) or not mortality_releases:
        raise ValueError('Declare mortality linkage release types')
    kinds=set(mortality_releases)
    if not kinds<= {'public_lmf2019','restricted_lmf2019','restricted_lmf2022'}:
        raise ValueError('Unreviewed mortality linkage release')
    if 'public_lmf2019' in kinds and len(kinds)>1:
        raise ValueError('Do not combine public and restricted mortality linkages')
    if len(kinds)>1:raise ValueError('Do not mix linkage vintages as one outcome dataset')
    return {'single_linkage_release':True,'data_access_authorized_by_this_check':False}


def precision_planning(events: int) -> dict:
    """Transparent event-count approximation, not EPV criterion or AUC inference.

    For an independent Poisson count, approximate relative SE=1/sqrt(D).
    Survey clustering/weights, censoring and fitted parameters are NOT included.
    Zero events cannot support this approximation or a separately fitted interval.
    """
    if type(events) is not int or events<0:raise ValueError('Nonnegative count required')
    return {'events':events,'poisson_relative_se_approx':1/math.sqrt(events) if events else None,
            'model_sample_size_criterion':False,'accounts_for_complex_survey':False}


def mec_weight_1999_2010(year: int, two_year: float, four_year: float | None=None) -> float:
    """CDC pooled-weight rule: the first4years use provided4year weights.

    This is a twelve-year examination-weight function, not the cardiac-protein
    subsample rule. Within single cycles use the original two-year weight.
    """
    if type(year) is not int or year not in range(1999,2010,2):
        raise ValueError('Expected one of six1999-2010 survey cycles')
    value=four_year if year<2003 else two_year
    if type(value) not in (int,float) or not math.isfinite(value) or value<=0:
        raise ValueError('Positive finite appropriate survey weight required')
    return value/3 if year<2003 else value/6


def profile_transfer_gate(required: Sequence[str], available: Sequence[str], *,
                          target: str, survey_start: int) -> dict:
    """Metadata-only gate; prevents pretending a cohort-wide absent assay is measured.

    Matching names do not establish assay equivalence or performance. Missing
    assays demand a newly fitted reduced model or an independently tested method.
    """
    for values in (required,available):
        if isinstance(values,(str,bytes)) or not isinstance(values,Sequence):
            raise ValueError('Feature names must be sequences')
        if any(not isinstance(v,str) or not v.strip() for v in values):
            raise ValueError('Nonempty feature names required')
        if len(set(values))!=len(values):raise ValueError('Duplicate feature definition')
    endpoint=observability(target,survey_start=survey_start)
    missing=sorted(set(required)-set(available))
    return {'same_model_replay_supported_by_metadata':not missing and endpoint['label_observable'],
            'missing_assays':missing,'endpoint':endpoint,'clinical_use_ready':False,
            'action':'fit_and_freeze_reduced_model_before_test_outcomes' if missing else
                     'obtain_sufficient_endpoint_detail' if not endpoint['label_observable'] else
                     'verify_measurement_methods_followup_and_independent_performance'}
