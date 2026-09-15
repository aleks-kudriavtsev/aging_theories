# Multicause20: scope before new data inspection

Date: 2026-09-15. Base commit: 1e5772ad813ef0d3f1a18f36c4bf34949972e039.

The user-supplied integrated RU/DE/US research plan defines twenty intended modules, not five: ischemic heart disease; heart failure/cardiomyopathy; cerebrovascular disease; hypertensive disease; atrial fibrillation/flutter; aortic aneurysm/dissection; chronic lower respiratory disease; Alzheimer/dementia; Parkinson disease; diabetes; renal disease; chronic liver disease/cirrhosis; lung; colorectal; pancreatic; breast; prostate; liver; gastric; hematologic malignancy mortality.

These are a target taxonomy, NOT twenty empirically identified outcomes in the existing recent NHANES public mortality files. The latter cannot be split by assigning national shares or treating missing subcategories as zero. Infection and external mortality are retained as competitors, not counted toward the target of at least ten noninfectious groups.

## Data-resolution strategy
Inspect original public documentation before participant outcomes. Older NHANES III UCOD_113 releases will only be used if a legitimate documented source is recovered. A second official candidate is NHEFS, whose 1992 public mortality file retains underlying ICD-9 causes. Its baseline biomarker and follow-up definitions must be reviewed before specifying a historical proof-of-method experiment. Restricted later linkage will not be accessed or reconstructed. Older US results cannot be presented as current RU/DE validation.

## Scientific invariants
- Train joint cause-specific intensities from observed individual measurements and mutually exclusive underlying causes. Sum CIFs with survival to one. Never multiply an all-cause score by national cause fractions.
- Keep target modules, source-observable endpoints and empirically estimated heads separately counted. Unknown death cause remains an event.
- Prespecify the cohort, split, feature blocks, penalization, horizons and primary loss in a second pre-performance analysis document after codebook inspection. No claim of preregistration for choices made after inspection.
- Expanded markers enter numerical fitting only if actually observed and harmonized. NT-proBNP is not BNP, hs-cTnT is not hs-cTnI, TNF-alpha is not sTNFR1. GDF15/IL6/GFAP/NT-proBNP and other repository nominees may be recorded as pending-measurement candidates but not assigned invented coefficients.
- No individual data, identifiers, institutional contracts, consent documents or proprietary aptamer sequences enter git. Public artifacts: code, aggregate results, source manifests, exact model definition and abstract.
- Report all prespecified comparisons. A positive claim in the abstract must be generated from executed results; no outcome-dependent favourable panel selection.

This is an auditable research decision log, not an ethics approval or clinical trial registration. Existing frozen models and the public application remain unchanged.
