"""Public CDC-only retrieval; resolves actual release-index hrefs before data."""
from __future__ import annotations
import argparse,concurrent.futures,datetime,json
from html.parser import HTMLParser
from pathlib import Path
import urllib.parse,urllib.request
from ..nfl22.analyze import sha
from .gate import require_qualification
ROOT=Path(__file__).parent
INDEX='https://wwwn.cdc.gov/nchs/nhanes/nhanes3/datafiles.aspx'
STATIC={'adult.sas':'https://wwwn.cdc.gov/nchs/data/nhanes3/1a/adult.sas',
        'cystatin.pdf':'https://wwwn.cdc.gov/nchs/data/nhanes3/27a/SSCYSTAT.pdf',
        'btp.htm':'https://wwwn.cdc.gov/nchs/data/nhanes3/36a/SSNH3BTP.htm',
        'mortality.dat':'https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_III_MORT_2019_PUBLIC.dat'}
DYNAMIC={'adult.dat':('/1a/','adult.dat'),'cystatin.xpt':('/27a/','sscystat.xpt'),'btp.xpt':('/36a/','ssnh3btp.xpt')}

class Links(HTMLParser):
    def __init__(self):super().__init__();self.hrefs=[]
    def handle_starttag(self,tag,attrs):
        if tag.lower()=='a':
            self.hrefs.extend(v for k,v in attrs if k.lower()=='href' and v)

def allowed_url(url):
    p=urllib.parse.urlsplit(url)
    return p.scheme=='https' and p.hostname in {'wwwn.cdc.gov','ftp.cdc.gov'} and p.username is None and p.password is None and p.port in (None,443)

def resolve_links(html):
    p=Links();p.feed(html);result={}
    for name,(directory,filename) in DYNAMIC.items():
        matches=set()
        for href in p.hrefs:
            url=urllib.parse.urljoin(INDEX,href)
            if urllib.parse.urlsplit(url).scheme.lower()=='ftp':
                parts=urllib.parse.urlsplit(url);url=urllib.parse.urlunsplit(('https',parts.netloc,parts.path,parts.query,parts.fragment))
            part=urllib.parse.urlsplit(url)
            if directory in part.path.lower() and part.path.lower().endswith('/'+filename) and allowed_url(url):matches.add(url)
        if len(matches)!=1:raise ValueError('Unresolved/ambiguous official release href: '+name)
        result[name]=matches.pop()
    return result

def retrieve(out):
    c=json.loads((ROOT/'CONTRACT.json').read_text());require_qualification(c)
    out=Path(out);out.mkdir(parents=True)
    def get(item):
        name,url=item
        if not allowed_url(url):raise ValueError('Nonofficial source')
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'MortalityCohort27/1.0'}),timeout=60) as response:
                if not allowed_url(response.url):raise ValueError('Unexpected resolved host')
                raw=response.read(100_000_001);final=response.url
            if len(raw)>100_000_000:raise ValueError('Oversize source')
            if name.endswith('.xpt') and not raw.startswith(b'HEADER RECORD'):raise ValueError('Not XPORT')
            if name.endswith('.dat') and b'<html' in raw[:1000].lower():raise ValueError('Not fixed-width data')
            if name.endswith('.pdf') and not raw.startswith(b'%PDF'):raise ValueError('Not PDF')
            path=out/name;path.write_bytes(raw)
            return dict(file=name,url=url,resolved_url=final,status='retrieved',bytes=len(raw),sha256=sha(path))
        except Exception as error:return dict(file=name,url=url,status='failed',error=str(error))
    first=get(('index.htm',INDEX));rows=[first]
    if first['status']!='retrieved':raise RuntimeError('Public metadata index unavailable')
    urls={**STATIC,**resolve_links((out/'index.htm').read_text(encoding='utf-8-sig'))}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:rows.extend(pool.map(get,urls.items()))
    doc=dict(plan_commit='a15532f',lock_sha256=sha(ROOT/'LOCK.json'),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=rows)
    (out/'manifest.json').write_text(json.dumps(doc,indent=2)+'\n')
    print(json.dumps(dict(retrieved=sum(x['status']=='retrieved' for x in rows),failures=[r for r in rows if r['status']!='retrieved'])))
    if any(r['status']!='retrieved' for r in rows):raise RuntimeError('Required official source unavailable')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    retrieve(p.parse_args().out)
