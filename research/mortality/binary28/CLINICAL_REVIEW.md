# NHANES III simple clinical core28 — primary metadata verification

2026-10-29. Metadata only: official CDC adult.sas and ADULT-acc.pdf were read.
No participant records, biomarker values, endpoints or new comparison results
were opened. No repository file changed.

## Selected core and interpretation

A separate NHANESIII-compatible eight-year binary model can use interview age,
sex, three-state smoking and the released overall systolic BP mean. No frozen
compact4/expanded13 history feature is being reconstructed. This core does not
adjust for creatinine/eGFR/UACR, full laboratory panels or every established
mortality risk factor; any protein gain is incremental to this limited core.

Ten selected numeric fields, one-based inclusive SAS INPUT positions:

| Field | Positions | Verified meaning / role |
|---|---|---|
| SEQN | 1–5 | Within-survey person linkage key; never exported |
| HSSEX | 15 | Reported sex code1 male /2 female |
| HSAGEIR | 18–19 | Screener/interview age quantity |
| HSAGEU | 20 | Age unit; require2(years), not months |
| SDPPSU6 | 43 | Combined six-year public pseudo-PSU1/2 |
| SDPSTRA6 | 44–45 | Combined six-year public pseudo-stratum01–49 |
| HAR1 | 2281 | Whether at least100 cigarettes over entire life |
| HAR3 | 2285 | Whether currently smokes cigarettes |
| HAZMNK1R | 3337–3339 | Overall arithmetic mean systolic/K1 BP, mmHg |
| HAZNOK1R | 3340–3341 | Number of measurements contributing to HAZMNK1R |

SAS positions were checked inside INPUT, not LENGTH or FORMAT. The metadata
uses separate variable storage lengths, which must not be interpreted as
source positions. Selected support now ends at3341, so every record must
contain at least3341characters; the qualification27 six-field minimum45 is
insufficient for this expanded reader. Preserve original bytes and every
record meeting the documented format; do not pad missing selected fields.

## Smoking: exact questionnaire categories and skips

Official codebook PDF page-index288 (displayed page289):
- HAR1:1=yes,2=no,8=blank but applicable.
- HAR3:1=yes,2=no,8=blank but applicable. It is structurally blank after
  HAR1=2, which skips to HAR14. This skip is not an unknown smoker state.

Prespecified three-state mapping:
- never:HAR1=2, regardless of structurally missing HAR3;
- current:HAR1=1 andHAR3=1;
- former:HAR1=1 andHAR3=2;
- any other HAR1/HAR3 combination:missing; do not default8/blank to never.

The question captures current cigarette smoking, not every nicotine/tobacco
product or pack-years. Do not relabel HAR3 as daily versus occasional status:
its released answer is simple yes/no.

## Systolic BP: exact released definition and sentinels

Official PDF page-index389: HAZMNK1R is overall average K1/systolic BP in mmHg;
reported range069–246;888 denotes blank but applicable and physical blank is
missing. HAZNOK1R has released codes01–06,88=blank but applicable, plus blank.
Require a known count1..6 for the selected BP feature;88/blank is not a large
count and does not establish a valid mean. The reported69–246 range is a
release-frequency range, not a universal clinical/instrument validity limit.
It can be asserted as source-release support; do not clip unknown values into it.

Official PDF page-index432 explicitly defines HAZMNK1R as the arithmetic mean
of up to six K1 readings: at most three from the initial household interview
and at most three from the initial MEC examination. If household measurements
are absent, initial MEC readings may supply the mean. The derived mean may be
odd-valued; do not round to an even number or recompute using a guessed rule.

This is a released summary of the FIRST survey interview/exam measurements.
It is not release37a's later repeated-examination convenience sample and does
not shift the biological/mortality baseline to the2006/2009 stored-serum assay
years. Retain original MEC-origin mortality followup. A household/MEC mixture
is an explicit measurement limitation and differs from a strictly MEC-only BP.

## Alternatives reviewed, excluded from this chosen simple core

Diabetes is available in adult without another source: HAD1col1561 is doctor-
diagnosed diabetes; borderline/prediabetic is coded2(no). Pregnancy distinction
requires HAD3col1563 and HAD4col1564 plus sex. A bare HAD1-positive variable
would include gestational history; do not silently call it non-gestational.
Parent selected core excludes diabetes to avoid this extra harmonization.

Measured BMI is in EXAM, not ADULT. Official exam.sas gives BMPBMI1524–1527,
FORMAT Z6.1; source also has measured weight/height and flags. Its nominal
record width6233(+CRLF), hence the full separate EXAM.dat is much larger than
ADULT; actual bytes would require metadata/retrieval inspection. No EXAM
participant file was fetched. Avoid estimating actual file size as verified.
No BMI or self-reported weight/height proxy is used in the selected core.

## Primary sources

- SAS layout/labels: https://wwwn.cdc.gov/nchs/data/nhanes3/1a/adult.sas
  INPUT selected fields around lines1390–1418,2252–2254,2624–2625;
  LABEL lines3493–3495 and3865–3866.
- Official codebook: https://wwwn.cdc.gov/nchs/data/nhanes3/1a/ADULT-acc.pdf
  PDF zero-based pages288(smoking),389(BP values/counts),392(age notes),
  394(design),432(BP derived-mean definition); pages170(diabetes).
- BMI alternative metadata only: https://wwwn.cdc.gov/nchs/data/nhanes3/1a/exam.sas

These metadata facts support a new clinical model specification, not clinical
validity or transport of any prior model's coefficients.
