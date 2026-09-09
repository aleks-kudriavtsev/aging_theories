"""Targeted XPORT zero repair from exact source bytes, NOT epsilon clipping.

The installed pandas XPORT decoder can decode IBM zero 0000000000000000 as
5.397605346934028e-79. Upstream report: pandas-dev/pandas#30051. We correct only
bitwise zero (including signed zero); every other decoded value is untouched.
Pandas still supplies the file layout and nonzero decoding. This is NOT a fully
independent SAS parser. No row-level or identifying data are emitted in audits.
"""
from __future__ import annotations
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd


def exact_zero_rows(raw: bytes, start: int, record_length: int, nrows: int,
                    position: int, length: int) -> np.ndarray:
    if (any(type(v) is not int for v in [start,record_length,nrows,position,length])
        or min(start,nrows,position)<0 or record_length<=0 or not 2<=length<=8
        or position+length>record_length or start+nrows*record_length>len(raw)):
        raise ValueError('Unsupported/truncated XPORT layout')
    a=np.ndarray(shape=(nrows,length),dtype='u1',buffer=raw,
                 offset=start+position,strides=(record_length,1))
    return np.all(a[:,1:]==0,axis=1)&((a[:,0]==0)|(a[:,0]==128))


def read_xpt_checked(path: Path) -> pd.DataFrame:
    path=Path(path);raw=path.read_bytes()
    if not raw.startswith(b'HEADER RECORD'):raise ValueError('Not an XPORT file')
    with pd.read_sas(path,format='xport',iterator=True) as reader:
        start,step=int(reader.record_start),int(reader.record_length)
        fields=reader.fields.copy();d=reader.read()
    changed={};zero_counts={}
    for field in fields:
        if field['ntype']!='numeric':continue
        name=field['name'].decode('ascii');position=int(field['npos']);length=int(field['field_length'])
        mask=exact_zero_rows(raw,start,step,len(d),position,length)
        zero_counts[name]=int(mask.sum())
        changed[name]=int((mask & (d[name].to_numpy()!=0)).sum())
        d.loc[mask,name]=0.
    if 'SEQN' not in d or d.SEQN.isna().any() or d.SEQN.duplicated().any():
        raise ValueError('Missing/nonunique participant key')
    d.attrs['xpt_zero_audit']={'filename':path.name,'sha256':hashlib.sha256(raw).hexdigest(),
      'nrows':len(d),'numeric_zero_counts':{k:v for k,v in zero_counts.items() if v},
      'corrected_counts':{k:v for k,v in changed.items() if v},'pandas_version':pd.__version__,
      'epsilon_threshold_used':False}
    return d
