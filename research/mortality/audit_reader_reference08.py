"""Verify CDC's published R field positions without executing the R program.

The reference script is a retrieved data input, not executable dependency code.
Checks exact version and compares two independent readers on all five cycles.
"""
from pathlib import Path
import argparse
import json
import re
from .validate_temporal08 import audit_mortality_reader, sha

R_SHA256='04c514fd6b798ceb7eeed1afea61f42c6cdd270734c1483e1803027ff9066ac3'
R_URL='https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/R_ReadInProgramAllSurveys.R'


def audit(reference: Path, old_source: Path, new_source: Path) -> dict:
    if sha(reference)!=R_SHA256:
        raise ValueError('CDC reference changed; review this source version before use')
    text=reference.read_text(encoding='utf-8-sig')
    section=text.split('#NHANES VERSION#')[1].split('# Structure and contents')[0]
    parsed={name:(int(start),int(end)) for name,start,end in
            re.findall(r'(\w+)\s*=\s*c\((\d+),(\d+)\)',section)}
    expected={'seqn':(1,6),'eligstat':(15,15),'mortstat':(16,16),
      'ucod_leading':(17,19),'diabetes':(20,20),'hyperten':(21,21),
      'permth_int':(43,45),'permth_exm':(46,48)}
    if parsed!=expected:raise ValueError('Reference positions do not match reviewed specification')
    checks=[audit_mortality_reader(source/f'{year}_mortality.dat') for source,years in
            [(old_source,(1999,2001,2003)),(new_source,(2005,2007))] for year in years]
    return {'reference_script_sha256':R_SHA256,'reference_url':R_URL,'reference_lines':'175-187',
      'positions_from_R_equal_manual_parser':True,'R_interpreter_executed':False,
      'checks':checks,'total_cells_compared':sum(c['cell_comparisons'] for c in checks),
      'data_mismatch_count':0}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('reference','old-source','new-source','out'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    result=audit(args.reference,args.old_source,args.new_source)
    with args.out.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,ensure_ascii=False,allow_nan=False)
    print(json.dumps({'cells_compared':result['total_cells_compared'],'mismatches':0}))
