# NHANES III BMI source review29 — 2026-10-29

Metadata only; no participant file was requested or opened.

## Definition and coding

The official first-examination codebook identifies BMPBMI as weight divided by squared height in metres. INPUT positions are1524–1527; released numeric BMI is in kg/m².8888 and blank mean missing. The documented7.3–79.6 span is a release frequency range, not a clinical exclusion rule. Read decimal text directly: Z6.1 is an OUTPUT FORMAT, not an implied-decimal INPUT instruction; do not divide arbitrary integer tokens by10.

Released weight/height may contain CDC substitutions. Preserve this qualification and preferably export aggregate source flags: BMPWTFLG1514 accepts0/3/4; BMPHTFLG1533 accepts0/1/2/3/5.0 means no substitution;3 denotes modelled reported information. This is not new study imputation. A predeclared all-released-BMI primary analysis avoids silently selecting a different cohort; primary BMI must not be labelled uniformly direct measured BMI. Clothing/refusal metadata does not supply a separate mandatory adult BMI exclusion rule here.

Source: https://wwwn.cdc.gov/nchs/data/nhanes3/1a/exam-acc.pdf , PDF indices197–205 (BMI/table/derivation and substitution notes).

## Safe reader and acquisition

https://wwwn.cdc.gov/nchs/data/nhanes3/1a/exam.sas separately contains LENGTH,FORMAT and INPUT blocks. Verify only INPUT positions; numeric storage length8 and output format are not source locations. FILENAME declares LRECL6235 including2CRLF positions; the nominal last INPUT endpoint is6233. INFILE uses MISSOVER. Do not require every physical line to retain unselected trailing spaces: selected support ends at1527, or1533 when flags are imported, with maximum nominal6233. Record actual width frequencies; do not skip arbitrary short records. Compare a manual decimal parser with independent pandas.read_fwf on identical selected positions.

The official codebook counts31,311 examination records (30,818MEC+493home; independently its BMI categories sum27,830+227+3,254). Maximum nominal CRLF payload is31,311×6,235=195,224,085bytes. A per-file220,000,000byte bound is metadata-justified and finite;100MB is too small. Prefer streamed download/hash rather than unbounded read. Check actual bytes/SHA and source count after acquisition; no actual data byte count was verified here.

## Exact public route

Official index https://wwwn.cdc.gov/nchs/nhanes/nhanes3/datafiles.aspx lists Examination File DAT190MB, SAS and the lowercase PDF URL above. Expected source basename is1a/exam.dat; resolve its actual href from this index before fetching. Local attempt to read ONLY index HTML was network403, so raw DAT href was not independently resolved. Never substitute3a/examse.dat or use the revised dietary files for BMI. Preserve WTCYPEX6 for the shared protein domain, not a general exam weight.
