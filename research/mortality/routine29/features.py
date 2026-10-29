"""Prespecified renal/routine features; source variables retain released values.

Calibration applies only to serum CEP. eGFR is the adult2021 creatinine-only
estimate, never a measured GFR or a cystatin-derived comparator. Published
missing/nondetect codes are interpreted before any formula or logarithm.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from ..binary28.cohort import CORE as SIMPLE_CORE

CORE = SIMPLE_CORE + ('bmi', 'hba1c', 'egfr_low', 'egfr_high', 'log_uacr')
ROUTINE = CORE


def _positive(values, missing_code):
    """Positive finite released numerical values, without frequency cutoffs."""
    return values.where(np.isfinite(values) & values.gt(0) & values.ne(missing_code))


def engineer(frame):
    """Return a copy with the fixed routine core and explicit provenance flags.

    ``frame`` has the six already engineered original clinical predictors plus
    HSAGEIR/HSSEX and raw CEP,GHP,GHPMETH,URP,UBP,BMPBMI,BMPWTFLG,BMPHTFLG.
    Nothing here changes or imputes the raw source columns. The caller has
    already restricted the survey/age/endpoint domain and confirmed units.
    """
    if not isinstance(frame, pd.DataFrame):
        raise ValueError('A source DataFrame is required')
    required = set(SIMPLE_CORE) | {'HSAGEIR', 'HSSEX', 'CEP', 'GHP', 'GHPMETH',
                                  'URP', 'UBP', 'BMPBMI', 'BMPWTFLG', 'BMPHTFLG'}
    if not required <= set(frame.columns):
        raise ValueError('Missing prespecified clinical/routine source fields')
    d = frame.copy()

    creatinine = _positive(d.CEP, 8888)
    calibrated = .960 * creatinine - .184
    d['creatinine_calibrated'] = calibrated.where(np.isfinite(calibrated) & calibrated.gt(0))
    female = d.HSSEX.eq(2)
    age = d.HSAGEIR.where(np.isfinite(d.HSAGEIR) & d.HSAGEIR.ge(18))
    ratio = d.creatinine_calibrated / np.where(female, .7, .9)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        egfr = (142 * np.minimum(ratio, 1) ** np.where(female, -.241, -.302)
                * np.maximum(ratio, 1) ** -1.200 * .9938 ** age
                * np.where(female, 1.012, 1.))
    d['egfr'] = egfr.where(d.HSSEX.isin([1, 2]) & np.isfinite(egfr) & egfr.gt(0))
    d['egfr_low'] = np.minimum(d.egfr, 60) / 15
    d['egfr_high'] = np.maximum(d.egfr - 60, 0) / 30

    urine_creatinine = _positive(d.URP, 88888)
    # CDC calls released7.9 statistically suspect/BLD; it is not a measured
    # denominator. The primary domain excludes all released URP below10.
    primary_urine_creatinine = urine_creatinine.where(urine_creatinine.ge(10))
    urine_albumin = _positive(d.UBP, 888888)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        ratio_uacr = 100 * urine_albumin / primary_urine_creatinine
    d['uacr'] = ratio_uacr.where(np.isfinite(ratio_uacr) & ratio_uacr.gt(0))
    d['log_uacr'] = np.log2(d.uacr)

    released_bmi = _positive(d.BMPBMI, 8888)
    d['bmi'] = (released_bmi - 27) / 5
    ghp = _positive(d.GHP, 8888)
    d['hba1c'] = ghp.where(d.GHPMETH.isin([1, 2, 3, 4]))

    # Flags are for aggregate reporting and do not add model predictors.
    d['urine_albumin_reported_below_detection'] = urine_albumin.eq(.4)
    d['urine_creatinine_reported_suspect'] = urine_creatinine.lt(10)
    d['serum_creatinine_source_missing'] = creatinine.isna()
    d['serum_creatinine_calibration_nonpositive'] = creatinine.notna() & d.creatinine_calibrated.isna()
    d['routine_missing_bmi'] = d.bmi.isna()
    d['routine_missing_hba1c'] = d.hba1c.isna()
    d['routine_missing_egfr'] = d.egfr.isna()
    d['routine_missing_uacr'] = d.uacr.isna()
    return d


def primary_complete(d):
    """Same11-finite-predictor mask for every compared routine-domain arm."""
    return np.isfinite(d[list(CORE)]).all(axis=1).to_numpy(dtype=bool)


def strict_metadata_mask(d):
    """Metadata sensitivity only; intersect with the PRIMARY common domain.

    No documented anthropometric substitution, primary HPLC, above the
    released albumin fill. This is not an individual QC/LOQ certification.
    """
    albumin = _positive(d.UBP, 888888)
    return (d.BMPWTFLG.eq(0) & d.BMPHTFLG.eq(0)
            & d.GHPMETH.isin([1, 2, 3]) & albumin.gt(.4)).to_numpy(dtype=bool)
