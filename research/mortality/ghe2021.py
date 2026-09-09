"""Read-only WHO GHE XLSX extraction and auditable, nonoverlapping rankings.

Standard-library ZIP/XML reader, not a general spreadsheet editor. The scientific
classification is an explicit, versioned policy. Raw NCD != strict noninfectious.
No individual risk, regression coefficient or clinical cutoff is inferred here.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, re, zipfile
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
SEX = {'Persons':'all', 'Males':'male', 'Females':'female'}
COUNTRIES = ('RUS','DEU','USA')
AGE_SHEETS = ['All ages','0-4','5-14','15-29','30-49','50-59','60-69','70+']
FILES = {
 'deaths': ('0_ghe2021_deaths_bycountry_2021.xlsx','1f593b791db8a9408fdce1973ef2b16b4cd5c9cf8e1838b661e0be9ffa5b7418'),
 'asdr': ('1_ghe2021_deathrates_bycountry_asdr.xlsx','3c87314a348384b7de86e4f5c3b5534637bb6f5982897df2812bd106c05cb863'),
 'cdr': ('2_ghe2021_deathrates_bycountry_cdr.xlsx','dbb22132411aa5b6191d51bcb88697ab6752416e1ecb50c75e62bbda7e090daa')}


def iter_rows(path: Path, sheet_index: int):
    """Read cached numeric/text values only. Formula cells are rejected."""
    with zipfile.ZipFile(path) as z:
        strings = [''.join(t.text or '' for t in el.iter(NS+'t'))
                   for el in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        with z.open(f'xl/worksheets/sheet{sheet_index}.xml') as stream:
            for _, row in ET.iterparse(stream, events=('end',)):
                if row.tag != NS+'row': continue
                result = {}
                for cell in row:
                    if cell.find(NS+'f') is not None:
                        raise ValueError('Formula cell in input; evaluate/review upstream')
                    ref = cell.get('r','')
                    col = re.sub(r'\d', '', ref)
                    v = cell.find(NS+'v')
                    if v is None: continue
                    value = strings[int(v.text)] if cell.get('t')=='s' else v.text
                    result[col] = value
                yield int(row.get('r')), result
                row.clear()


def _finite(value: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0: raise ValueError('Invalid nonnegative estimate')
    return result


def read_sheet(path: Path, sheet_index: int, *, metric: str, year: int, age: str):
    records, nodes, population = [], {}, {}
    columns, stack = {}, []
    for row, cells in iter_rows(path, sheet_index):
        if row == 8:
            columns = {v:k for k,v in cells.items() if v in COUNTRIES}
            if set(columns) != set(COUNTRIES): raise ValueError('Country header changed')
        if cells.get('A') not in SEX: continue
        sex = SEX[cells['A']]
        if 'B' not in cells:
            if 'Population' in cells.get('D',''):
                for iso, col in columns.items():
                    population[(iso,sex)] = 1000*_finite(cells[col])
            continue
        code = int(cells['B'])
        label_col = max(c for c in ('D','E','F','G') if c in cells)
        name = cells[label_col].replace('\xa0',' ')
        depth = ord(label_col)-ord('D')
        if code == 0: depth = -1
        # Excel indents 1700 under injuries; official Annex A identifies IV,
        # and all-cause reconciliation requires it to be a separate root child.
        if code == 1700: depth = 0
        if sex == 'all':
            while stack and nodes[stack[-1]]['depth'] >= depth: stack.pop()
            parent = stack[-1] if stack else None
            nodes[code] = {'code':code,'name_en':name,'depth':depth,
                           'parent':parent,'children':[]}
            if parent is not None: nodes[parent]['children'].append(code)
            stack.append(code)
        for iso,col in columns.items():
            value = _finite(cells[col])
            if metric == 'deaths': value *= 1000  # Source explicitly in thousands.
            records.append({'country':iso,'sex':sex,'age':age,'year':year,
                            'ghe_code':code,'metric':metric,'value':value,
                            'population':population[(iso,sex)],
                            'source_file':path.name,'source_sheet_index':sheet_index,
                            'source_cell':f'{col}{row}'})
    return records, nodes


def descendants(nodes: dict, code: int) -> set[int]:
    result = set()
    for child in nodes[code]['children']:
        result.add(child); result.update(descendants(nodes,child))
    return result


def disease_cut(nodes: dict) -> list[int]:
    """Complete partition at disease-group depth; retain shallower leaves."""
    cut = [code for code,node in nodes.items() if code and
           (node['depth']==2 or (not node['children'] and node['depth']<2))]
    assert_antichain(nodes,cut)
    return cut


def assert_antichain(nodes: dict, codes: list[int]):
    if len(codes)!=len(set(codes)): raise ValueError('Duplicate cause')
    for code in codes:
        if code not in nodes: raise ValueError('Unknown cause')
        if descendants(nodes,code).intersection(codes):
            raise ValueError('Parent and child would be double counted')


def check_close(actual: float, expected: float, description: str) -> float:
    error = abs(actual-expected)
    if not math.isclose(actual,expected,rel_tol=1e-8,abs_tol=1e-6):
        raise ValueError(f'{description}: {actual} != {expected}')
    return error


def audit(records, nodes, cut):
    """Fail rather than silently publish unreconciled counts/rates."""
    index = {}; stats = defaultdict(int); maxerr = defaultdict(float)
    for row in records:
        key = tuple(row[k] for k in ('country','sex','age','year','ghe_code','metric'))
        if key in index: raise ValueError('Duplicate observation')
        index[key] = row
    strata = {key[:4]+(key[5],) for key in index}
    for iso,sex,age,year,metric in strata:
        def val(code): return index[(iso,sex,age,year,code,metric)]['value']
        err = check_close(val(0),math.fsum(val(c) for c in cut),'complete partition')
        stats['complete_partition_checks']+=1; maxerr['partition']=max(maxerr['partition'],err)
        for code,node in nodes.items():
            if node['children']:
                err=check_close(val(code),math.fsum(val(c) for c in node['children']),'hierarchy')
                stats['hierarchy_checks']+=1; maxerr['hierarchy']=max(maxerr['hierarchy'],err)
    for iso in COUNTRIES:
        for age in AGE_SHEETS:
            for code in nodes:
                v=lambda sex:index[(iso,sex,age,2021,code,'deaths')]['value']
                err=check_close(v('all'),v('male')+v('female'),'sex sum')
                stats['sex_sum_checks']+=1;maxerr['sex_sum']=max(maxerr['sex_sum'],err)
        for sex in SEX.values():
            for code in nodes:
                val=lambda age:index[(iso,sex,age,2021,code,'deaths')]['value']
                err=check_close(val('All ages'),math.fsum(val(a) for a in AGE_SHEETS[1:]),'age sum')
                stats['age_sum_checks']+=1;maxerr['age_sum']=max(maxerr['age_sum'],err)
                row=index[(iso,sex,'All ages',2021,code,'deaths')]
                expected=row['value']/row['population']*100000
                actual=index[(iso,sex,'All ages',2021,code,'cdr')]['value']
                err=check_close(actual,expected,'CDR and count/population')
                stats['cdr_checks']+=1;maxerr['cdr']=max(maxerr['cdr'],err)
    return {'status':'passed','observations':len(records),'tree_nodes':len(nodes),
            'partition_groups':len(cut),'checks':dict(stats),'max_absolute_errors':dict(maxerr),
            'note':'Arithmetic/data-integrity audit, not external clinical validation'}, index


def rank_bounds(counts: dict[int,float], compatible: list[int], mixed: list[int]):
    """Conditional ranking bounds, not confidence intervals.

    Known compatible counts fixed; each mixed group may contribute 0..its total.
    Fixed grouping; no redistribution of pandemic residual or WHO uncertainty.
    Equal counts use ascending GHE code as a deterministic tie-breaker.
    """
    if len(compatible)!=len(set(compatible)) or len(mixed)!=len(set(mixed)):
        raise ValueError('Duplicate policy cause')
    for code in compatible+mixed:
        if code not in counts or counts[code]<0 or not math.isfinite(counts[code]):
            raise ValueError('Invalid count')
    if set(compatible)&set(mixed): raise ValueError('Overlapping policy classes')
    ordered=sorted(compatible,key=lambda c:(-counts[c],c))
    out=[]
    for rank,code in enumerate(ordered,1):
        worse=sum((counts[m]>counts[code] or (counts[m]==counts[code] and m<code)) for m in mixed)
        out.append({'ghe_code':code,'rank_compatible':rank,'rank_lower':rank,'rank_upper':rank+worse})
    return out


def write_csv(path: Path, rows: list[dict]):
    if not rows: raise ValueError('Cannot export empty table')
    with path.open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def build(source_dir: Path, policy_path: Path, out: Path):
    out.mkdir(parents=True,exist_ok=True)
    policy=json.loads(policy_path.read_text(encoding='utf-8'))
    records=[]; nodes=None
    for metric,(filename,sha) in FILES.items():
        path=source_dir/filename
        if hashlib.sha256(path.read_bytes()).hexdigest()!=sha:
            raise ValueError(f'Source version changed: {filename}; review instead of silent reuse')
        work=[(i+2,2021,age)for i,age in enumerate(AGE_SHEETS)] if metric=='deaths' else [(2,2021,'All ages'),(4,2019,'All ages')]
        for sheet,year,age in work:
            new,newnodes=read_sheet(path,sheet,metric=metric,year=year,age=age)
            if nodes is None: nodes=newnodes
            if newnodes!=nodes: raise ValueError('Cause hierarchy mismatch between sheets')
            records.extend(new)
    cut=disease_cut(nodes)
    if set(map(int,policy['causes']))!=set(cut):raise ValueError('Policy does not cover complete partition')
    review,index=audit(records,nodes,cut)
    compatible=[c for c in cut if policy['causes'][str(c)]['eligibility']=='compatible']
    mixed=[c for c in cut if policy['causes'][str(c)]['eligibility']=='mixed_or_unreviewed']
    ncd=descendants(nodes,600).intersection(cut)
    rankings=[];raw=[];summary=[];sex_age=[];sensitivity=[]
    for iso in COUNTRIES:
        def rec(code,metric='deaths',sex='all',age='All ages',year=2021):
            return index[(iso,sex,age,year,code,metric)]
        counts={c:rec(c)['value'] for c in cut}
        bounds=rank_bounds(counts,compatible,mixed)
        for r in bounds[:20]:
            c=r['ghe_code'];p=policy['causes'][str(c)];d=rec(c)
            rankings.append({'country':iso,**r,'name_ru':p['name_ru'],'name_en':nodes[c]['name_en'],
                'deaths_estimate':d['value'],'share_all_deaths':d['value']/rec(0)['value'],
                'cdr_per_100000':rec(c,'cdr')['value'],'asdr_per_100000':rec(c,'asdr')['value'],
                'icd10':p.get('icd10'),'year':2021,'source_sheet':'All ages','source_cell':d['source_cell']})
            for sex in SEX.values():
                for age in AGE_SHEETS:
                    d=rec(c,sex=sex,age=age)
                    sex_age.append({'country':iso,'ghe_code':c,'sex':sex,'age':age,'year':2021,
                        'deaths_estimate':d['value'],'population_estimate':d['population'],
                        'age_specific_or_crude_rate_per_100000':d['value']/d['population']*100000,
                        'source_cell':d['source_cell'],'source_sheet':age})
            for year in (2019,2021):
                for sex in SEX.values():
                    sensitivity.append({'country':iso,'ghe_code':c,'year':year,'sex':sex,
                         'asdr_per_100000':rec(c,'asdr',sex=sex,year=year)['value'],
                         'cdr_per_100000':rec(c,'cdr',sex=sex,year=year)['value']})
        for rank,c in enumerate(sorted(ncd,key=lambda c:(-counts[c],c))[:20],1):
            p=policy['causes'][str(c)]
            raw.append({'country':iso,'rank_raw_ncd':rank,'ghe_code':c,'name_ru':p['name_ru'],
                        'name_en':nodes[c]['name_en'],'deaths_estimate':counts[c],
                        'asdr_per_100000':rec(c,'asdr')['value'],'eligibility':p['eligibility'],
                        'reason':p['reason']})
        byclass={cl:math.fsum(counts[c]for c in cut if policy['causes'][str(c)]['eligibility']==cl)
                 for cl in set(p['eligibility']for p in policy['causes'].values())}
        summary.append({'country':iso,'all_deaths':rec(0)['value'],'population':rec(0)['population'],
                        'raw_ncd_deaths':rec(600)['value'],**byclass,
                        'compatible_top20_sum':sum(r['deaths_estimate']for r in rankings if r['country']==iso),
                        'other_pandemic_unallocated':rec(1700)['value']})
    write_csv(out/'top20_compatible.csv',rankings);write_csv(out/'top20_raw_ghe_ncd.csv',raw)
    write_csv(out/'top20_sex_age.csv',sex_age);write_csv(out/'top20_rates_2019_2021.csv',sensitivity)
    write_csv(out/'country_audit.csv',summary)
    write_csv(out/'source_excerpt.csv',records)
    (out/'data_audit.json').write_text(json.dumps(review,ensure_ascii=False,indent=2))
    (out/'rankings.json').write_text(json.dumps({'status':'research_conditional_group_rankings_not_individual_risk',
        'source_version':'WHO GHE2021 release2024; files retrieved2026-09-09','year':2021,
        'top20_compatible':rankings,'top20_raw_ncd':raw,'countries':summary,'audit':review},ensure_ascii=False,indent=2))
    return review

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources',type=Path,required=True)
    parser.add_argument('--policy',type=Path,default=Path(__file__).with_name('ghe2021_policy.json'))
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args(); print(json.dumps(build(args.sources,args.policy,args.out),indent=2))
