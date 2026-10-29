# Representation amendment, before any qualification endpoint counts

2026-10-29, after failed run40887539474, before first successful person-count/
eight-year endpoint analysis. ANALYSIS_PLAN.md and scientific inclusion rules
remain unchanged. This records a parser correction, not a new scientific sample.

Official adult.sas contains separate LENGTH and INPUT blocks. Numeric storage
lengths such as SEQN7 are not byte positions SEQN1–5. The reader now isolates
INPUT and compares complete numeric ranges; regression fixtures cover LENGTH
collisions, CRLF, and spaced malformed ranges (HSSEX15–16).

The failed real reader produced ONLY aggregate source-width diagnostics:
3346characters:20014records;3345characters:36records, total20050. The official
SAS program uses INFILE ADULT ... MISSOVER, which permits absent trailing
fields.3346 is the last nominal INPUT column, not a guarantee of every physical
line's length. All six imported fields lie within columns1–45.

Accordingly, when the trusted source SAS declares MISSOVER, the selected-field
reader accepts physical lengths45…3346, preserving the original source bytes
and EVERY record; without that declaration exact3346 is still required.
No record is padded/truncated or silently skipped; shorter-than45/overlong
records, wrong INPUT positions, nonnumeric selected fields or duplicate IDs
still stop. An independent pandas FWF reader must exactly match all120300
selected cells, and the official adult release n20050 is checked before joins.
The source-width histogram and MISSOVER evidence are preserved in results.

No outcomes, lab concentrations, comparative metrics or participant records
were printed by diagnostics. Tests cover accepted trailing absence, rejected
selected-field truncation and the absence of a MISSOVER declaration.

Official sources:
https://wwwn.cdc.gov/nchs/data/nhanes3/1a/adult.sas
https://wwwn.cdc.gov/nchs/data/nhanes3/1a/ADULT-acc.pdf
The metadata-only independent review is NHANESIII_SOURCE_AUDIT.md.
All NCHS surveys share LMF layout; examination followup for NHANESIII starts
at MEC date: https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf (pp1,3).
