"""Audit the licensed public Olink file without fabricating mortality records.

Dube, Raw Olink data, doi:10.6084/m9.figshare.22811915, CC BY 4.0.
NPX stays NPX. Repeated panels remain separate OlinkIDs. Only aggregates leave
this procedure; raw public participant identifiers are not written to outputs.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import re
import zipfile

ZIP_SHA256='602c708713529c5e85f41640a8781b46e35abc44e3e48e7282ee71ffbb855ac6'
ZIP_MD5='3fa86946ab0099a3be44dca8871c0e29'
TARGETS={'GDF15':'Q99988','IL6':'P05231','TNF':'P01375','GFAP':'P14136',
         'NEFL':'P07196','CST3':'P01034','NPPB_unresolved_proteoform':'P16860'}
FOUR=('GDF15','IL6','TNF','GFAP')
PATTERN=re.compile(r'^SSNA-([0-9]+[A-Z]?)-(PR1|PT1|PR2|PT2)$')
REQUIRED={'SampleID','OlinkID','UniProt','Assay','Panel','Panel_Lot_Nr','PlateID',
          'QC_Warning','LOD','NPX','Normalization','Assay_Warning'}


def parse_sample(sample):
    match=PATTERN.fullmatch(sample)
    return match.groups() if match else None


def quantitative_qc(row):
    """A conservative usable-above-LOD flag, NOT an independently verified LOQ."""
    try:npx,lod=float(row['NPX']),float(row['LOD'])
    except (KeyError,TypeError,ValueError):return {'numeric':False,'below_LOD':None,'usable_above_LOD':False}
    numeric=math.isfinite(npx) and math.isfinite(lod)
    below=npx<lod if numeric else None
    return {'numeric':numeric,'below_LOD':below,'usable_above_LOD':numeric and not below
            and row.get('QC_Warning')=='PASS' and row.get('Assay_Warning')=='PASS'}


def audit_rows(rows):
    subjects=set();samples=set();visits={};controls=set();seen=set();assays={};all_assays=set()
    nrows=0;biological_rows=0
    uniprot_to_id={v:k for k,v in TARGETS.items()}
    for row in rows:
        nrows+=1
        if not REQUIRED<=row.keys():raise ValueError('Missing mandatory Olink metadata')
        pair=(row['SampleID'],row['OlinkID'])
        if pair in seen:raise ValueError('Duplicate sample/OlinkID; do not average replicates silently')
        seen.add(pair);all_assays.add(row['OlinkID'])
        sample=parse_sample(row['SampleID'])
        if sample is None:controls.add(row['SampleID']);continue
        subject,visit=sample;subjects.add(subject);samples.add(row['SampleID']);visits.setdefault(subject,set()).add(visit)
        biological_rows+=1
        target=uniprot_to_id.get(row['UniProt'])
        if target is None:continue
        identity=(target,row['OlinkID'],row['UniProt'],row['Assay'],row['Panel'],row['Panel_Lot_Nr'],row['Normalization'])
        code=row['OlinkID']
        if code not in assays:
            assays[code]={'identity':identity,'n':0,'below':0,'nonnumeric':0,'qc_failed':0,'assay_warning':0,
                          'usable':set(),'baseline_usable':set(),'plate_ids':set()}
        a=assays[code]
        if a['identity']!=identity:raise ValueError('An OlinkID has inconsistent analytical metadata')
        a['n']+=1;a['plate_ids'].add(row['PlateID']);q=quantitative_qc(row)
        a['nonnumeric']+=not q['numeric'];a['below']+=q['below_LOD'] is True
        a['qc_failed']+=row['QC_Warning']!='PASS';a['assay_warning']+=row['Assay_Warning']!='PASS'
        if q['usable_above_LOD']:
            a['usable'].add(row['SampleID'])
            if visit=='PR1':a['baseline_usable'].add(row['SampleID'])
    table=[]
    for code,a in sorted(assays.items()):
        target,oid,uniprot,label,panel,lot,normalization=a['identity']
        table.append({'target':target,'OlinkID':oid,'UniProt':uniprot,'source_Assay':label,'panel':panel,
          'lot':lot,'normalization':normalization,'unit_scale':'NPX','material':'not_confirmed_in_download',
          'sample_timepoints':a['n'],'below_LOD':a['below'],'not_finite':a['nonnumeric'],
          'sample_QC_warning':a['qc_failed'],'assay_warning':a['assay_warning'],
          'PASS_and_at_or_above_LOD':len(a['usable']),'baseline_PR1_PASS_above_LOD':len(a['baseline_usable']),
          'plate_ids':'|'.join(sorted(a['plate_ids'])),
          'exact_NT_proBNP_confirmed':False if target=='NPPB_unresolved_proteoform' else None})
    options={t:[oid for oid,a in assays.items() if a['identity'][0]==t] for t in FOUR}
    intersections=[]
    if all(options.values()):
        for ids in itertools.product(*(options[t] for t in FOUR)):
            good=set.intersection(*(assays[oid]['usable'] for oid in ids))
            baseline=set.intersection(*(assays[oid]['baseline_usable'] for oid in ids))
            intersections.append({**dict(zip(FOUR,ids)),'joint_PASS_above_LOD_timepoints':len(good),
              'joint_PASS_above_LOD_subjects':len({parse_sample(x)[0] for x in good}),
              'joint_baseline_PR1_PASS_above_LOD':len(baseline),
              'not_selected_by_outcome':True})
    summary={'raw_rows':nrows,'biological_rows':biological_rows,'assay_ids':len(all_assays),
      'biological_subjects':len(subjects),'biological_timepoints':len(samples),'control_sample_ids':len(controls),
      'every_subject_has_all_four_visits':all(v=={'PR1','PT1','PR2','PT2'} for v in visits.values()),
      'requested_four_measured':all(options.values()),'cross_panel_assays_not_averaged':True,
      'joint_four_assay_combinations':len(intersections),
      'max_joint_four_PASS_above_LOD_timepoints':max((r['joint_PASS_above_LOD_timepoints'] for r in intersections),default=0),
      'max_joint_four_baseline_PASS_above_LOD':max((r['joint_baseline_PR1_PASS_above_LOD'] for r in intersections),default=0),
      'mortality_followup_in_files':False,'expanded13_clinical_profile_in_files':False,
      'mortality_training_ready':False,'NT_PROBNP_proteoform_confirmed':False,
      'concentration_mass_units':False,'clinical_use_ready':False,
      'participant_ids_exported':False,'criterion':'NPX finite, NPX >= LOD, both warnings PASS; not a validated LOQ'}
    return summary,table,intersections


def run(archive,out):
    archive,out=Path(archive),Path(out)
    if out.exists():raise ValueError('New output directory required')
    blob=archive.read_bytes()
    if hashlib.sha256(blob).hexdigest()!=ZIP_SHA256 or hashlib.md5(blob).hexdigest()!=ZIP_MD5:
        raise ValueError('Reviewed Figshare source version changed')
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        members=[m for m in z.infolist() if m.filename.endswith('.csv')]
        if len(members)!=2 or sum(m.file_size for m in members)>50_000_000:raise ValueError('Unexpected archive layout')
        fingerprints=[];rows=[]
        for member in sorted(members,key=lambda m:m.filename):
            raw=z.read(member)
            fingerprints.append({'file':Path(member.filename).name,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
            reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')),delimiter=';')
            rows.extend(reader)
    summary,table,intersections=audit_rows(rows)
    summary.update(source_doi='10.6084/m9.figshare.22811915',source_file_id=40553531,license='CC BY 4.0',
        author='Marie-Pierre Dube',zip_sha256=ZIP_SHA256,source_files=fingerprints,
        experimental_context='10 participants, heat acclimation, four repeated timepoints')
    out.mkdir(parents=True)
    for name,data in [('target_assays.csv',table),('four_protein_joint_qc.csv',intersections)]:
        with (out/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    (out/'JOINT_DATA_AUDIT.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2));return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.archive,a.out)
