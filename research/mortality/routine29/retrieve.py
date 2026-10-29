"""Public NHANES III original LAB/EXAM extension with metadata-bounded sizes."""
from __future__ import annotations
import argparse,concurrent.futures,datetime,json
from pathlib import Path
import urllib.parse,urllib.request
from ..cohort27 import retrieve as prior
from ..nfl22.analyze import sha
ROOT=Path(__file__).parent
STATIC={'lab.sas':'https://wwwn.cdc.gov/nchs/data/nhanes3/1a/lab.sas',
        'exam.sas':'https://wwwn.cdc.gov/nchs/data/nhanes3/1a/exam.sas',
        'lab.pdf':'https://wwwn.cdc.gov/nchs/data/nhanes3/1a/lab-acc.pdf',
        'exam.pdf':'https://wwwn.cdc.gov/nchs/data/nhanes3/1a/exam-acc.pdf'}
DYNAMIC={'lab.dat':'lab.dat','exam.dat':'exam.dat'}


def check_study():
    lock=json.loads((ROOT/'LOCK.json').read_text())
    for name,digest in lock['study_metadata_sha256'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Prespecified routine29 metadata changed')
    if sha(ROOT.parent/'binary28/SOURCE_MANIFEST.json')!=lock['prior28_source_manifest_sha256']:
        raise ValueError('Prior28 source manifest changed')
    return lock


def resolve_additional(html):
    parser=prior.Links();parser.feed(html);result={}
    for name,filename in DYNAMIC.items():
        matches=set()
        for href in parser.hrefs:
            url=urllib.parse.urljoin(prior.INDEX,href);parts=urllib.parse.urlsplit(url)
            if parts.scheme.lower()=='ftp':url=urllib.parse.urlunsplit(('https',parts.netloc,parts.path,parts.query,parts.fragment))
            parts=urllib.parse.urlsplit(url)
            if '/1a/' in parts.path.lower() and parts.path.lower().endswith('/'+filename) and prior.allowed_url(url):matches.add(url)
        if len(matches)!=1:raise ValueError('Unresolved/ambiguous original-exam source href: '+name)
        result[name]=matches.pop()
    return result


def retrieve(out):
    check_study();prior.retrieve(out) # old public permission and source gates retained
    out=Path(out);base=json.loads((out/'manifest.json').read_text())
    urls={**STATIC,**resolve_additional((out/'index.htm').read_text(encoding='utf-8-sig'))}
    def get(item):
        name,url=item;cap=220_000_000 if name=='exam.dat' else 100_000_000
        if not prior.allowed_url(url):raise ValueError('Nonofficial source')
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'MortalityRoutine29/1.0'}),timeout=90) as response:
                if not prior.allowed_url(response.url):raise ValueError('Unexpected resolved host')
                data=response.read(cap+1);final=response.url
            if len(data)>cap:raise ValueError('Metadata-bounded source size exceeded')
            if name.endswith('.pdf') and not data.startswith(b'%PDF'):raise ValueError('Not PDF')
            if not name.endswith('.pdf') and b'<html' in data[:1000].lower():raise ValueError('Not declared SAS/text data')
            path=out/name;path.write_bytes(data)
            return dict(file=name,url=url,resolved_url=final,status='retrieved',bytes=len(data),sha256=sha(path),maximum_download_bytes=cap)
        except Exception as error:return dict(file=name,url=url,status='failed',error=str(error))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:extra=list(pool.map(get,urls.items()))
    doc=dict(plan_commit='7800458',lock_sha256=sha(ROOT/'LOCK.json'),base_manifest_sha256=sha(out/'manifest.json'),
             retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=base['files']+extra)
    (out/'routine29_manifest.json').write_text(json.dumps(doc,indent=2)+'\n')
    print(json.dumps(dict(additional_retrieved=sum(x['status']=='retrieved' for x in extra),failures=[x for x in extra if x['status']!='retrieved'])))
    if any(x['status']!='retrieved' for x in extra):raise RuntimeError('New official LAB/EXAM source unavailable')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    retrieve(p.parse_args().out)
