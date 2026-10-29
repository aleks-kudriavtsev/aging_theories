"""Retrieve public CDC inputs; fail if scientific sources are unavailable."""
import urllib.request,urllib.parse,pathlib,hashlib,json,concurrent.futures,datetime
import argparse
p=argparse.ArgumentParser();p.add_argument('--out',type=pathlib.Path,required=True)
root=p.parse_args().out;root.mkdir(parents=True)
components=['DEMO','BPX','BMX','SMQ','MCQ','BPQ','DIQ','BIOPRO','TCHOL','HDL','GHB','ALB_CR','CBC','SSSNFL']
urls={f'2013_{c}_H.{ext}':f'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/{c}_H.{ext}' for c in components for ext in ['xpt','htm']}
urls['2013_mortality.dat']='https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_2013_2014_MORT_2019_PUBLIC.dat'
def get(item):
    name,url=item
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'NfL22-Research/1.0'}),timeout=50) as r:
            if urllib.parse.urlparse(r.url).hostname not in {'wwwn.cdc.gov','ftp.cdc.gov'}:raise ValueError('Unexpected host')
            data=r.read(20_000_001)
        if len(data)>20_000_000:raise ValueError('Size exceeded')
        if name.endswith('.xpt') and not data.startswith(b'HEADER RECORD'):raise ValueError('Not XPORT')
        if name.endswith('.dat') and b'<html' in data[:1000].lower():raise ValueError('Not mortality')
        (root/name).write_bytes(data)
        return {'file':name,'url':url,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data),'status':'retrieved'}
    except Exception as exc:return {'file':name,'url':url,'status':'failed','error':str(exc)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(get,urls.items()))
(root/'manifest.json').write_text(json.dumps({'source':'CDC public only','date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':rows},indent=2))
print(json.dumps({'retrieved':sum(r['status']=='retrieved' for r in rows),'failures':[r for r in rows if r['status']!='retrieved']},indent=2))
if any(r['status']!='retrieved' and r['file'].endswith(('.dat','.xpt')) for r in rows):raise SystemExit('Required source unavailable')