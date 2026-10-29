"""Public-only source intake for a protocol-fixed exploratory B2M comparison."""
from __future__ import annotations
import argparse
import concurrent.futures
import datetime
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import urllib.parse
import urllib.request

ROOT = Path(__file__).parent
PLAN_COMMIT = '6774fbe'
HOSTS = {'wwwn.cdc.gov', 'www.cdc.gov', 'ftp.cdc.gov'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_study():
    lock = json.loads((ROOT / 'LOCK.json').read_text())
    for name, expected in lock['metadata_sha256'].items():
        path = ROOT.parent / name if '/' in name else ROOT / name
        if sha(path) != expected:
            raise ValueError('Locked study metadata changed')
    for name, expected in lock['historical_scientific_python_sha256'].items():
        if sha(ROOT.parent / name) != expected:
            raise ValueError('Historical scientific code changed')
    return lock


def allowed_url(url):
    try:
        p = urllib.parse.urlsplit(url)
        return (p.scheme == 'https' and p.hostname in HOSTS and p.username is None
                and p.password is None and p.port in (None, 443) and bool(p.path)
                and not p.fragment)
    except (ValueError, TypeError):
        return False


class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == 'a' and dict(attrs).get('href'):
            self.hrefs.append(dict(attrs)['href'])


def resolve_sscard(html, base):
    parser = Links(); parser.feed(html); matches = set()
    for href in parser.hrefs:
        url = urllib.parse.urljoin(base, href)
        if (allowed_url(url) and Path(urllib.parse.urlsplit(url).path).name.lower() == 'sscard_a.xpt'):
            matches.add(url)
    if len(matches) != 1:
        raise ValueError('Unique official SSCARD_A data anchor unresolved; no guessed URL fallback')
    return matches.pop()


def static_urls():
    spec = json.loads((ROOT / 'SOURCE_SPEC.json').read_text()); urls = {}
    for cycle in spec['cycles']:
        yr = cycle['year']
        for component, row in cycle['components'].items():
            urls[f'{yr}_{component}.xpt'] = row['historical_declared_data_url']
            urls[f'{yr}_{component}.htm'] = row['documentation_url']
        urls[f'{yr}_mortality.dat'] = cycle['mortality_url']
    urls.update({'1999_SSCARD_A.htm': spec['sscard']['metadata_url'],
        'weighting.htm': spec['weighting_url'], 'variance.htm': spec['variance_url'],
        'mortality_dictionary.pdf': spec['mortality_dictionary_url'],
        'mortality_description.pdf': spec['mortality_description_url']})
    return urls


def get_one(out, name, url):
    if not allowed_url(url):
        raise ValueError('Nonofficial or unsafe source URL')
    with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'MortalityB2M31/1.0'}), timeout=90) as r:
        if not allowed_url(r.url):
            raise ValueError('Nonofficial resolved URL')
        data = r.read(25_000_001); final = r.url
    if len(data) > 25_000_000:
        raise ValueError('Source size bound exceeded')
    if name.endswith('.xpt') and not data.startswith(b'HEADER RECORD'):
        raise ValueError('Not declared XPORT data')
    if name.endswith('.pdf') and not data.startswith(b'%PDF'):
        raise ValueError('Not declared PDF documentation')
    if name.endswith('.htm') and b'<html' not in data[:10000].lower():
        raise ValueError('Not declared HTML metadata')
    if name.endswith('.dat') and b'<html' in data[:10000].lower():
        raise ValueError('Not declared mortality text')
    (out / name).write_bytes(data)
    return dict(file=name, url=url, resolved_url=final, status='retrieved', bytes=len(data), sha256=sha(out/name))


def retrieve(out):
    check_study(); out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    spec = json.loads((ROOT/'SOURCE_SPEC.json').read_text()); rows = []
    try:
        index = get_one(out, 'index.htm', spec['sscard_index_url']); rows.append(index)
        sscard = resolve_sscard((out/'index.htm').read_text(encoding='utf-8-sig'), index['resolved_url'])
        urls = {**static_urls(), '1999_SSCARD_A.xpt': sscard}
        def get(item):
            name, url = item
            try:
                return get_one(out, name, url)
            except Exception as error:
                return dict(file=name, url=url, status='failed', error=str(error))
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            rows.extend(pool.map(get, urls.items()))
    except Exception as error:
        rows.append(dict(file='source_resolution', status='failed', error=str(error)))
    doc = dict(schema='b2m31-retrieval/1', plan_commit=PLAN_COMMIT, lock_sha256=sha(ROOT/'LOCK.json'),
               retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), files=rows)
    (out/'manifest.json').write_text(json.dumps(doc, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(retrieved=sum(r['status']=='retrieved' for r in rows),
                         failures=[r for r in rows if r['status']!='retrieved'])))
    if any(r['status'] != 'retrieved' for r in rows):
        raise RuntimeError('Official source intake blocked before comparison')
    return doc


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--out', type=Path, required=True)
    retrieve(p.parse_args().out)
