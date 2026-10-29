"""New11-term core from Continuous NHANES1999–2004 released measurements.

This is a source-specific research core, not the frozen NHANESIII adapter.
``RIDAGEYR`` is screening age. The caller supplies ``raw_creatinine`` in mg/dL
from current LAB18.LBXSCR, L40_B.LBDSCR or L40_C.LBXSCR respectively and checks
the source columns/units before engineering. CDC's current1999 analytic note
requires1.013*original+.147;2001–04 require no additional standardization.
Source: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/LAB18.htm

SBP here is the prespecified research mean of positive finite BPXSY1..3, not an
official report-average variable. Urine quantities remain reported values:
albumin ug/mL and creatinine mg/dL. No NHANESIII serum correction, urine floor,
HbA1c regression, cystatin-derived eGFR, clipping or imputation is applied.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


CORE = ('age_linear', 'age_quadratic', 'male', 'current', 'former', 'sbp',
        'bmi', 'hba1c', 'egfr_low', 'egfr_high', 'log_uacr')
SUPPORTED_CYCLES = (1999, 2001, 2003)
REQUIRED = ('cycle', 'RIDAGEYR', 'RIAGENDR', 'SMQ020', 'SMQ040',
            'BPXSY1', 'BPXSY2', 'BPXSY3', 'BMXBMI', 'LBXGH',
            'raw_creatinine', 'URXUMA', 'URXUCR')


def _numeric(values):
    """Retain finite numeric source values; booleans are not measurements."""
    numeric = pd.to_numeric(values, errors='coerce').astype(float)
    boolean = values.map(lambda value: isinstance(value, (bool, np.bool_)))
    return numeric.where(np.isfinite(numeric) & ~boolean)


def _positive(values):
    numeric = _numeric(values)
    return numeric.where(numeric.gt(0))


def engineer(frame):
    """Copy input and append fixed predictors plus source-validity flags.

    Unsupported/missing cycles fail because no serum mapping is approved for
    them. Other missing or invalid source values yield missing dependent terms,
    so ``primary_complete`` defines the common domain without imputation. The
    caller still qualifies participants, endpoint, design and assay provenance.
    No person-level export or fitting occurs here.
    """
    if not isinstance(frame, pd.DataFrame):
        raise ValueError('A source DataFrame is required')
    if not set(REQUIRED) <= set(frame.columns):
        raise ValueError('Missing prespecified continuous-NHANES source fields')
    d = frame.copy()
    cycle = _numeric(d.cycle)
    if not cycle.isin(SUPPORTED_CYCLES).all():
        raise ValueError('Unsupported or missing serum-creatinine cycle mapping')

    screening_age = _numeric(d.RIDAGEYR)
    valid_age = (screening_age.between(40, 75)
                 & screening_age.eq(np.floor(screening_age)))
    age = screening_age.where(valid_age)
    d['age_linear'] = (age - 60) / 10
    d['age_quadratic'] = d.age_linear ** 2
    sex = _numeric(d.RIAGENDR)
    d['male'] = sex.map({1.: 1., 2.: 0.})

    smoked100 = _numeric(d.SMQ020)
    smoking_now = _numeric(d.SMQ040)
    never = smoked100.eq(2)
    current = smoked100.eq(1) & smoking_now.isin([1, 2])
    former = smoked100.eq(1) & smoking_now.eq(3)
    known_smoking = never | current | former
    d['current'] = np.where(known_smoking, current.astype(float), np.nan)
    d['former'] = np.where(known_smoking, former.astype(float), np.nan)

    systolic = pd.concat([_positive(d[name]) for name in ('BPXSY1', 'BPXSY2', 'BPXSY3')], axis=1)
    d['sbp_measurement_count'] = systolic.notna().sum(axis=1).astype(int)
    d['sbp'] = (systolic.mean(axis=1) - 130) / 20
    d['bmi'] = (_positive(d.BMXBMI) - 27) / 5
    # LBXGH is the original released DCCT-aligned HPLC percent; no regression.
    d['hba1c'] = _positive(d.LBXGH)

    serum = _positive(d.raw_creatinine)
    with np.errstate(over='ignore', invalid='ignore'):
        calibrated = serum.where(cycle.ne(1999), 1.013 * serum + .147)
    d['creatinine_calibrated'] = calibrated.where(np.isfinite(calibrated) & calibrated.gt(0))
    female = sex.eq(2)
    ratio = d.creatinine_calibrated / np.where(female, .7, .9)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        egfr = (142 * np.minimum(ratio, 1) ** np.where(female, -.241, -.302)
                * np.maximum(ratio, 1) ** -1.200 * .9938 ** age
                * np.where(female, 1.012, 1.))
    d['egfr'] = egfr.where(sex.isin([1, 2]) & np.isfinite(egfr) & egfr.gt(0))
    d['egfr_low'] = np.minimum(d.egfr, 60) / 15
    d['egfr_high'] = np.maximum(d.egfr - 60, 0) / 30

    urine_albumin = _positive(d.URXUMA)
    urine_creatinine = _positive(d.URXUCR)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        ratio_uacr = 100 * urine_albumin / urine_creatinine
    d['uacr'] = ratio_uacr.where(np.isfinite(ratio_uacr) & ratio_uacr.gt(0))
    d['log_uacr'] = np.log2(d.uacr)

    # Reportable source diagnostics only; these flags are not extra predictors.
    d['source_invalid_screening_age'] = ~valid_age
    d['source_invalid_sex'] = ~sex.isin([1, 2])
    d['smoking_unknown'] = ~known_smoking
    d['serum_creatinine_source_missing'] = serum.isna()
    d['serum_creatinine_calibration_applied'] = cycle.eq(1999)
    d['serum_creatinine_calibration_nonpositive'] = serum.notna() & d.creatinine_calibrated.isna()
    d['routine_missing_sbp'] = d.sbp.isna()
    d['routine_missing_bmi'] = d.bmi.isna()
    d['routine_missing_hba1c'] = d.hba1c.isna()
    d['routine_missing_egfr'] = d.egfr.isna()
    d['routine_missing_uacr'] = d.uacr.isna()
    return d


def primary_complete(frame):
    """Common finite11-term domain for every compared study arm."""
    return np.isfinite(frame[list(CORE)]).all(axis=1).to_numpy(dtype=bool)
