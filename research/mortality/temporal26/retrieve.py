"""Retrieve only public official CDC2017–2018 inputs into a new directory."""
from __future__ import annotations
import argparse
import concurrent.futures
import datetime
import json
from pathlib import Path
import urllib.request
import urllib.parse
from ..nfl22.analyze import COMPONENTS, sha

ROOT = Path(__file__).parent

def source_urls():
    urls = {f'2017_{c}_J.{ext}': f'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/{c}_J.{ext}'
            for c in COMPONENTS for ext in ('xpt', 'htm')}
    urls['2017_mortality.dat'] = 'https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_2017_2018_MORT_2019_PUBLIC.dat'
    return urls


def retrieve(out):
    out = Path(out); out.mkdir(parents=True)
    def get(item):
        name, url = item
        try:
            req = urllib.request.Request(url, headers={'User-Agent':'MortalityTemporal26/1.0'})
            with urllib.request.urlopen(req, timeout=50) as response:
                final = response.url
                if urllib.parse.urlparse(final).hostname not in {'wwwn.cdc.gov','ftp.cdc.gov'}:
                    raise ValueError('Unexpected resolved source host')
                data = response.read(20_000_001)
            if len(data)>20_000_000: raise ValueError('Oversize source')
            if name.endswith('.xpt') and not data.startswith(b'HEADER RECORD'): raise ValueError('Not XPORT')
            if name.endswith('.dat') and b'<html' in data[:1000].lower(): raise ValueError('Not mortality text')
            if name.endswith('.htm') and b'<html' not in data[:3000].lower(): raise ValueError('Not documentation')
            path=out/name; path.write_bytes(data)
            return dict(file=name,url=url,resolved_url=final,status='retrieved',bytes=len(data),sha256=sha(path))
        except Exception as error:
            return dict(file=name,url=url,status='failed',error=str(error))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        rows=list(pool.map(get,source_urls().items()))
    doc=dict(analysis_lock_commit='5fbb68c',lock_sha256=sha(ROOT/'LOCK.json'),
             retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=rows)
    (out/'manifest.json').write_text(json.dumps(doc,indent=2)+'\n')
    print(json.dumps({'retrieved':sum(r['status']=='retrieved' for r in rows),'failures':[r for r in rows if r['status']!='retrieved']}))
    if any(r['status']!='retrieved' for r in rows): raise RuntimeError('Required official source unavailable')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--out',type=Path,required=True)
    retrieve(p.parse_args().out)
