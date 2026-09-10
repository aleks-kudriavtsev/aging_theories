"""Registered mortality harmonization; not an individual risk model.

Source: WHO Mortality Database, 23 February 2026, documentation table 8.
A common tabulation is necessary because Russia reports list 101 while DE/US
report 104. Remainders use explicit parent-complement semantics, not a claim
that every printed legacy ICD range is exhaustive under later revisions.
Only national, common-year records enter the comparison. No downloaded code.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, re, zipfile
from collections import Counter, defaultdict
from pathlib import Path

VERSION='0.10'
EXPECTED_ARCHIVES={'mort_documentation71f9e29d-7e3f-41e6-aafc-c4c1775c7aa3.zip': 'c4a1db856a08621fd4cc6a8ee986707fa37181b0730d16da13964f41d80368e2', 'mort_availability.zip': '78e5004f1f4b393d96a7d661d3b6859d431e02f498c6cad21d92216df97a4205', 'mort_country_codes.zip': '8c410820356fc572845b5281b36f638e044a565f808c4e72efc8fb69b07df6b2', 'mort_notes.zip': 'b9f9fdad0e70d3b12d775411248a56abafec999c77a0b782c57022dcad806e84', 'mort_pop.zip': '948a71a91c0d877e9b762982f21a12519d54b12cbcb46f6a18eae7f51868952d', 'morticd10_part5.zip': '6baa52d694e95b1bc5c97c2905ea29e3502e43cc7e6bee8610e4c10862ac0d6f', 'morticd10_part6.zip': 'cfbb3374a460c71dfb943bef8556396335f5449eccf655adbb4575c5ee6b7197'}
COUNTRIES={'4272':('RUS','Россия'),'4085':('DEU','Германия'),'2450':('USA','США')}
STATES={'compatible','mixed_or_unresolved','infection','external','ill_defined'}
AGES=['0-4']+[f'{a}-{a+4}' for a in range(5,85,5)]+['85+','unknown']
# code, label, operational eligibility, scope. No code is inferred from counts.
DEFS=[
 ('1001','Инфекционные и паразитарные болезни','infection','A00-B99'),
 ('1027','ЗНО губы, полости рта и глотки','compatible','C00-C14'),
 ('1028','ЗНО пищевода','compatible','C15'),('1029','ЗНО желудка','compatible','C16'),
 ('1030','ЗНО ободочной, прямой кишки и ануса','compatible','C18-C21'),
 ('1031','ЗНО печени и внутрипечёночных желчных протоков','compatible','C22'),
 ('1032','ЗНО поджелудочной железы','compatible','C25'),('1033','ЗНО гортани','compatible','C32'),
 ('1034','ЗНО трахеи, бронхов и лёгкого','compatible','C33-C34'),
 ('1035','Меланома кожи','compatible','C43'),('1036','ЗНО молочной железы','compatible','C50'),
 ('1037','ЗНО шейки матки','compatible','C53'),('1038','ЗНО других и неуточнённых частей матки','compatible','C54-C55'),
 ('1039','ЗНО яичника','compatible','C56'),('1040','ЗНО предстательной железы','compatible','C61'),
 ('1041','ЗНО мочевого пузыря','compatible','C67'),('1042','ЗНО головного мозга, оболочек и других частей ЦНС','compatible','C70-C72'),
 ('1043','Неходжкинские лимфомы: историческая группа','compatible','C82-C85'),
 ('1044','Множественная миелома и плазмоклеточные ЗНО','compatible','C90'),
 ('1045','Лейкозы','compatible','C91-C95'),
 ('1046','Прочие злокачественные новообразования','compatible','C00-C97 minus explicitly listed tumour groups'),
 ('1047','Новообразования, кроме злокачественных','compatible','D00-D48'),
 ('1049','Анемии','compatible','D50-D64'),('1050','Прочие болезни крови и иммунные нарушения','mixed_or_unresolved','D50-D89 minus D50-D64'),
 ('1052','Сахарный диабет','compatible','E10-E14'),('1053','Недостаточность питания','compatible','E40-E46'),
 ('1054','Прочие эндокринные, алиментарные и метаболические болезни','compatible','E00-E99 minus E10-E14,E40-E46'),
 ('1056','Психические расстройства вследствие употребления психоактивных веществ','compatible','F10-F19'),
 ('1057','Прочие психические и поведенческие расстройства','compatible','F00-F99 minus F10-F19'),
 ('1059','Менингит: выбранные коды','mixed_or_unresolved','G00,G03'),
 ('1060','Болезнь Альцгеймера','compatible','G30'),
 ('1061','Прочие болезни нервной системы','mixed_or_unresolved','G00-G99 minus G00,G03,G30'),
 ('1062','Болезни глаза и придаточного аппарата','mixed_or_unresolved','H00-H59'),
 ('1063','Болезни уха и сосцевидного отростка','mixed_or_unresolved','H60-H99'),
 ('1065','Острая ревматическая лихорадка и ревматические болезни сердца','mixed_or_unresolved','I00-I09'),
 ('1066','Гипертензивные болезни: I10-I13','compatible','I10-I13'),
 ('1067','Ишемическая болезнь сердца','compatible','I20-I25'),
 ('1068','Другие болезни сердца','mixed_or_unresolved','I26-I51'),
 ('1069','Цереброваскулярные болезни','compatible','I60-I69'),
 ('1070','Атеросклероз','compatible','I70'),
 ('1071','Прочие болезни системы кровообращения','mixed_or_unresolved','I00-I99 minus explicit circulatory groups'),
 ('1073','Грипп: историческая группа J10-J11','infection','J10-J11'),
 ('1074','Пневмония','infection','J12-J18'),('1075','Другие острые инфекции нижних дыхательных путей','infection','J20-J22'),
 ('1076','Хронические болезни нижних дыхательных путей','compatible','J40-J47'),
 ('1077','Прочие болезни органов дыхания','mixed_or_unresolved','J00-J99 minus explicit respiratory groups'),
 ('1079','Язвенная болезнь желудка и двенадцатиперстной кишки','compatible','K25-K27'),
 ('1080','Болезни печени','mixed_or_unresolved','K70-K76'),
 ('1081','Прочие болезни органов пищеварения','mixed_or_unresolved','K00-K99 minus K25-K27,K70-K76'),
 ('1082','Болезни кожи и подкожной клетчатки','mixed_or_unresolved','L00-L99'),
 ('1083','Болезни костно-мышечной системы и соединительной ткани','mixed_or_unresolved','M00-M99'),
 ('1085','Гломерулярные и тубулоинтерстициальные болезни почек','mixed_or_unresolved','N00-N15'),
 ('1086','Прочие болезни мочеполовой системы','mixed_or_unresolved','N00-N99 minus N00-N15'),
 ('1088','Беременность с абортивным исходом','mixed_or_unresolved','O00-O07'),
 ('1089','Другие прямые акушерские причины','mixed_or_unresolved','O10-O92'),
 ('1090','Непрямые акушерские причины','mixed_or_unresolved','O98-O99'),
 ('1091','Прочие материнские причины','mixed_or_unresolved','O00-O99 minus explicit maternal groups'),
 ('1092','Отдельные перинатальные состояния','mixed_or_unresolved','P00-P99'),
 ('1093','Врождённые аномалии и хромосомные нарушения','compatible','Q00-Q99'),
 ('1094','Симптомы и неточно определённые причины','ill_defined','R00-R99'),
 ('1095','Внешние причины','external','V00-Y99'),
 ('1901','Тяжёлый острый респираторный синдром','infection','U04.9'),
 ('1902','Нарушение, связанное с вейпингом','mixed_or_unresolved','U07.0'),
 ('1903','COVID-19','infection','U07.1-U07.2'),
 ('other_special','Прочие специальные или неотображённые коды','mixed_or_unresolved','unmapped special codes')]
POLICY={c:{'code':c,'label_ru':l,'status':s,'definition':d,'residual':('Прочие' in l or 'Другие болезни' in l)} for c,l,s,d in DEFS}
PARENTS={'1026':[str(c) for c in range(1027,1048)],'1048':['1049','1050'],
 '1051':['1052','1053','1054'],'1055':['1056','1057'],'1058':['1059','1060','1061'],
 '1064':[str(c) for c in range(1065,1072)],'1072':[str(c) for c in range(1073,1078)],
 '1078':['1079','1080','1081'],'1084':['1085','1086'],'1087':['1088','1089','1090','1091']}
# Complements close the partition; do not erase source's retained original codes.
RANGES={
 **{c:d for c,l,s,d in DEFS if re.fullmatch(r'[A-Z]\d\d(?:-[A-Z]\d\d)?',d)},
 '1059':'G00,G03'}
RANGES.pop('1001',None);RANGES.pop('1095',None)

def numeric_code(code):
    if not isinstance(code,str) or not re.fullmatch(r'[A-Z]\d{2}',code):raise ValueError('Expected three-character ICD-10 code')
    return (ord(code[0])-65)*100+int(code[1:])

def in_range(code,expression):
    value=numeric_code(code)
    for part in expression.split(','):
        lo,_,hi=part.partition('-')
        if numeric_code(lo)<=value<=numeric_code(hi or lo):return True
    return False

def map_detailed(code):
    """Exact assignment for ICD-10 103/104; no substrings of condensed codes."""
    if not isinstance(code,str) or not re.fullmatch(r'[A-Z]\d{2}[A-Z0-9]?',code):
        raise ValueError('Invalid detailed ICD-10: '+str(code))
    cat=code[:3];letter=cat[0]
    if letter in 'AB':return '1001'
    if letter in 'VWXY':return '1095'
    if code in {'U049'}:return '1901'
    if code=='U070':return '1902'
    if code in {'U071','U072'}:return '1903'
    # More specific leaves before chapter remainders.
    explicit=['1027','1028','1029','1030','1031','1032','1033','1034','1035','1036','1037','1038','1039','1040','1041','1042','1043','1044','1045','1047','1049','1052','1053','1056','1059','1060','1065','1066','1067','1068','1069','1070','1073','1074','1075','1076','1079','1080','1085','1088','1089','1090']
    for key in explicit:
        if in_range(cat,RANGES[key]):return key
    if letter=='C':return '1046'
    if letter=='D':return '1050'
    if letter=='E':return '1054'
    if letter=='F':return '1057'
    if letter=='G':return '1061'
    if letter=='H':return '1062' if int(cat[1:])<60 else '1063'
    return {'I':'1071','J':'1077','K':'1081','L':'1082','M':'1083','N':'1086','O':'1091','P':'1092','Q':'1093','R':'1094'}.get(letter,'other_special')

def integer(value):
    if not isinstance(value,str) or not re.fullmatch(r'\d+',value):raise ValueError('Nonnegative integer required, missing is not zero: '+repr(value))
    return int(value)

def age_vector(row,prefix='Deaths'):
    """0-4, 5-year bands through80-84,85+,unknown. Infant sub-bands not added."""
    fmt=str(row['Frmat']).zfill(2)
    if fmt not in {'00','01'}:raise ValueError('Unsupported age format '+fmt)
    cols=list(range(2,24))+(list(range(24,26)) if fmt=='00' else [])+[26]
    # In official ASCII, empty cells within an active age field denote zero;
    # require the complete age sum to equal the separately present total.
    values={i:integer(row.get(prefix+str(i)) or '0') for i in cols}
    if fmt=='01' and any(row.get(prefix+str(i),'') not in {'','0'} for i in [24,25]):
        raise ValueError('85+ format has additional occupied older-age cells')
    vec=[sum(values[i] for i in range(2,7))]+[values[i] for i in range(7,23)]+[sum(values.get(i,0) for i in range(23,26)),values[26]]
    total=integer(row[prefix+'1'])
    if sum(vec)!=total and prefix!='Pop':raise ValueError('Age counts do not reconcile with all-age total')
    if prefix=='Pop' and abs(sum(vec)-total)>12:raise ValueError('Population age totals differ beyond declared rounding tolerance')
    return vec

def add(a,b):return [x+y for x,y in zip(a,b)]

def national(row):return row.get('Admin1','')=='' and row.get('SubDiv','')==''

def validate_selection(rows):
    if not rows:raise ValueError('Empty country-year selection')
    keys={(r['Country'],r['Year'],r['List']) for r in rows}
    if len(keys)!=1:raise ValueError('Different country/year/tabulation in one selection')
    if not all(national(r) for r in rows):raise ValueError('Regional records not allowed in national comparison')
    listing=rows[0]['List']
    if listing not in {'101','103','104','10M'}:raise ValueError('Unreviewed tabulation')
    seen=set()
    for r in rows:
        if r['Sex'] not in {'1','2','9'}:raise ValueError('Unsupported sex coding')
        k=(r['Sex'],r['Cause'])
        if k in seen:raise ValueError('Duplicate cause-sex row')
        seen.add(k);age_vector(r)
    if listing in {'103','104','10M'}:
        for sex in set(r['Sex'] for r in rows):
            codes={r['Cause'] for r in rows if r['Sex']==sex and r['Cause']!='AAA'}
            if any(len(c)==4 and c[:3] in codes for c in codes):
                raise ValueError('Parent and child detailed ICD codes overlap')
    return listing

def harmonize(rows):
    listing=validate_selection(rows);out={};audit=[]
    totals={r['Sex']:age_vector(r) for r in rows if r['Cause']==('1000' if listing=='101' else 'AAA')}
    if not {'1','2'}.issubset(totals):raise ValueError('Missing sex-specific all-cause totals')
    for sex in totals:
        subset=[r for r in rows if r['Sex']==sex];raw={r['Cause']:r for r in subset}
        result={c:[0]*len(AGES) for c in POLICY}
        if listing=='101':
            unknown=set(raw)-set(POLICY)-set(PARENTS)-{'1000'}-{str(i) for i in range(1002,1026)}-{str(i) for i in range(1096,1104)}
            if unknown:raise ValueError('Unrecognized condensed codes: '+repr(unknown))
            for key in POLICY:
                if key in raw:result[key]=age_vector(raw[key])
            # Absence is accepted as zero only if parent AND all-cause arithmetic
            # close exactly; provenance reports omitted zero cells explicitly.
            for parent,children in PARENTS.items():
                expected=age_vector(raw[parent]) if parent in raw else [0]*len(AGES)
                actual=[sum(result[c][a] for c in children) for a in range(len(AGES))]
                if actual!=expected:raise ValueError('Condensed hierarchy does not reconcile: '+parent)
            for parent,children in [('1001',range(1002,1026)),('1095',range(1096,1104))]:
                expected=age_vector(raw[parent]) if parent in raw else [0]*len(AGES)
                actual=[sum(age_vector(raw[str(c)])[a] for c in children if str(c) in raw) for a in range(len(AGES))]
                if actual!=expected:raise ValueError('Excluded branch hierarchy does not reconcile: '+parent)
            audit.append({'sex':sex,'omitted_leaf_rows':sorted(set(POLICY)-set(raw)),
                          'omitted_zero_inference':'nonnegative counts and exact parent/all-cause identities'})
        else:
            for r in subset:
                if r['Cause']=='AAA':continue
                key=map_detailed(r['Cause']);result[key]=add(result[key],age_vector(r))
            audit.append({'sex':sex,'external_aggregate_only':len([r for r in subset if r['Cause'].startswith(tuple('VWXY'))])==1})
        actual=[sum(v[a] for v in result.values()) for a in range(len(AGES))]
        if actual!=totals[sex]:raise ValueError('Partition does not equal all causes by age and sex')
        out[sex]=result
    return out,totals,audit

def ranks_and_bounds(counts):
    """Conditional competition ranks at this fixed grouping; not confidence intervals.

    One unresolved group can contribute at most one rankable compatible subtotal.
    No redistribution of R-codes or splitting into more diseases is assumed.
    """
    certain={k:v for k,v in counts.items() if POLICY[k]['status']=='compatible'}
    uncertain={k:v for k,v in counts.items() if POLICY[k]['status'] in {'mixed_or_unresolved','ill_defined'}}
    return {k:{'rank_identified':1+sum(o>v for o in certain.values()),
               'rank_upper_at_fixed_grouping':1+sum(o>v for o in certain.values())+sum(o>v for o in uncertain.values())}
            for k,v in certain.items()}

def endpoint_relation(left,right):
    """Code-defined relation, never a transfer of coefficients to a subset."""
    if not left or not right:return 'unresolved'
    a,b=set(left),set(right)
    if a==b:return 'equal_code_set'
    if a<b:return 'proper_subset_no_coefficient_transfer'
    if a>b:return 'proper_superset_no_coefficient_transfer'
    return 'overlap_not_equivalent' if a&b else 'disjoint'

def code_set(*ranges):
    return {f'{chr(letter)}{n:02d}' for letter in range(65,91) for n in range(100)
            if any(in_range(f'{chr(letter)}{n:02d}',r) for r in ranges)}

def relation_audit():
    return [
      {'requested':'IHD','available':'NHANES_UCOD001_heart','relation':endpoint_relation(code_set('I20-I25'),code_set('I00-I09','I11','I13','I20-I51'))},
      {'requested':'all_dementia','available':'WHO_List1_Alzheimer','relation':endpoint_relation(code_set('F01-F03','G30'),code_set('G30'))},
      {'requested':'lung_cancer','available':'NHANES_UCOD002_all_malignant','relation':endpoint_relation(code_set('C33-C34'),code_set('C00-C97'))},
      {'requested':'COPD','available':'WHO_List1_chronic_lower_respiratory','relation':endpoint_relation(code_set('J40-J44'),code_set('J40-J47'))},
      {'requested':'IHD','available':'WHO_List1_1067','relation':endpoint_relation(code_set('I20-I25'),code_set('I20-I25'))}]

def read_sources(source):
    manifest=json.loads((source/'manifest.json').read_text(encoding='utf-8'))
    rows=[];scanned=0;verified=[]
    for item in manifest['files']:
        if item['status']!='retrieved':
            if item['file'] in EXPECTED_ARCHIVES:raise ValueError('Required reviewed archive not retrieved')
            continue
        name=item['file']
        if Path(name).name!=name:raise ValueError('Unsafe source filename')
        path=source/name;digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if digest!=item['sha256'] or digest!=EXPECTED_ARCHIVES.get(name) or path.stat().st_size!=item['bytes']:raise ValueError('Source version mismatch '+name)
        verified.append(name)
        if name not in {'morticd10_part5.zip','morticd10_part6.zip'}:continue
        with zipfile.ZipFile(path) as z:
            names=[n for n in z.namelist() if not n.endswith('/')]
            if len(names)!=1:raise ValueError('Unexpected source archive structure')
            with z.open(names[0]) as f:
                for line,r in enumerate(csv.DictReader(io.TextIOWrapper(f,encoding='utf-8-sig')),2):
                    scanned+=1
                    if r['Country'] in COUNTRIES and national(r):
                        rows.append({**r,'source_file':name,'source_line':line})
    if set(verified)!=set(EXPECTED_ARCHIVES):raise ValueError('Required standard source not verified')
    return rows,{'rows_scanned':scanned,'national_rows_selected':len(rows),'archives_verified':verified,'source_manifest_sha256':hashlib.sha256((source/'manifest.json').read_bytes()).hexdigest()}

def populations(source,year):
    result={}
    with zipfile.ZipFile(source/'mort_pop.zip') as z:
        with z.open(z.namelist()[0]) as f:
            for r in csv.DictReader(io.TextIOWrapper(f,encoding='utf-8-sig')):
                if r['Country'] in COUNTRIES and r['Year']==str(year) and national(r):
                    k=(r['Country'],r['Sex'])
                    if k in result:raise ValueError('Multiple population denominator rows')
                    ages=age_vector(r,'Pop')
                    result[k]={'ages':ages,'total':integer(r['Pop1']),'age_sum_minus_reported_total':sum(ages)-integer(r['Pop1'])}
    return result

def dump(path,rows):
    if not rows:return
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def build(source:Path,out:Path):
    if out.exists():raise ValueError('New output directory required')
    raw,audit=read_sources(source)
    available={c:sorted({int(r['Year']) for r in raw if r['Country']==c}) for c in COUNTRIES}
    shared=set.intersection(*(set(v) for v in available.values()))
    if not shared:raise ValueError('No common observed year')
    year=max(shared);pop=populations(source,year)
    ledger=[];age_rows=[];summary=[];raw_scope=[];reconciliation=[];detail=[]
    for country,(iso,label) in COUNTRIES.items():
        rows=[r for r in raw if r['Country']==country and int(r['Year'])==year]
        groups,totals,checks=harmonize(rows)
        sexgroups={**groups,'all':{k:[sum(groups[s][k][i] for s in groups) for i in range(len(AGES))] for k in POLICY}}
        totaltotal={**totals,'all':[sum(totals[s][i] for s in totals) for i in range(len(AGES))]}
        reconciliation.extend({'country':iso,**x} for x in checks)
        for sex,causes in sexgroups.items():
            total=sum(totaltotal[sex]);counts={k:sum(v) for k,v in causes.items()};ranks=ranks_and_bounds(counts)
            entry=pop.get((country,sex));den=entry['ages'] if entry else None
            reported_pop=entry['total'] if entry else None
            if sex=='all' and all((country,s) in pop for s in totals):
                den=[sum(pop[country,s]['ages'][i] for s in totals) for i in range(len(AGES))]
                reported_pop=sum(pop[country,s]['total'] for s in totals)
            # No US denominator is imported from an unrelated year. One checked
            # national Census denominator is reported only for total/all-age.
            all_pop=reported_pop if den is not None else 328239523 if (iso=='USA' and year==2019 and sex=='all') else None
            pop_source='WHO_MDB_same_year_population' if den is not None else 'US_Census_vintage2019_July1' if all_pop is not None else 'unavailable_same_year_stratum'
            for key,values in causes.items():
                base={'country':iso,'country_ru':label,'year':year,'sex':sex,'code':key,**{k:v for k,v in POLICY[key].items() if k!='code'},'deaths':sum(values),'share_all_deaths':sum(values)/total,'population':all_pop,'crude_rate_per_100k':sum(values)/all_pop*1e5 if all_pop else None,'population_source':pop_source,'source_tabulation':rows[0]['List'],**ranks.get(key,{'rank_identified':None,'rank_upper_at_fixed_grouping':None})}
                ledger.append(base)
                for i,age in enumerate(AGES):
                    n=den[i] if den is not None and age!='unknown' else None
                    age_rows.append({'country':iso,'year':year,'sex':sex,'code':key,'age':age,'deaths':values[i],'population':n,'age_specific_rate_per_100k':values[i]/n*1e5 if n else None,'status':POLICY[key]['status']})
            state=Counter()
            for k,v in counts.items():state[POLICY[k]['status']]+=v
            summary.append({'country':iso,'year':year,'sex':sex,'all_deaths':total,**{s:state[s] for s in sorted(STATES)},'unknown_age_deaths':totaltotal[sex][-1],'population':all_pop,'crude_all_rate_per_100k':total/all_pop*1e5 if all_pop else None,'population_source':pop_source,'ihd_plus_cerebrovascular_share':(counts['1067']+counts['1069'])/total})
        raw_scope.append({'country':iso,'source_tabulation':rows[0]['List'],'common_year':year,'observed_standard_years':available[country],'count_rows_in_common_year':len(rows),'latest_standard_observed_year':max(available[country]),'global_latest_year_not_established':True})
        # Granularity audit: retain selected detailed components, never backfill RU.
        for name,expr in [('dementia_F01_F03','F01-F03'),('Alzheimer_G30','G30'),('CKD_N18','N18'),('acute_kidney_N17','N17'),('infective_endocarditis_I33','I33'),('acute_myocarditis_I40','I40'),('liver_abscess_K750','K750')]:
            value=None
            detail_status='not_identifiable_from_condensed_tabulation'
            if rows[0]['List']=='101' and expr=='G30':
                value=sum(integer(r['Deaths1']) for r in rows if r['Cause']=='1060')
                detail_status='exact_condensed_group'
            if rows[0]['List']!='101':
                value=sum(integer(r['Deaths1']) for r in rows if r['Cause']!='AAA' and ((r['Cause']=='K750') if expr=='K750' else in_range(r['Cause'][:3],expr)))
            detail.append({'country':iso,'year':year,'component':name,'definition':expr,'deaths':value,'status':('direct_code_count' if rows[0]['List']!='101' else detail_status)})
    # All national counts are reconstructed, not adopted as test answers.
    primary_control={'RUS':1798307,'DEU':939520,'USA':2854838}
    actual={r['country']:r['all_deaths'] for r in summary if r['sex']=='all'}
    if year==2019 and actual!=primary_control:raise ValueError('Frozen extraction control differs; inspect version')
    out.mkdir(parents=True)
    dump(out/'harmonized_ledger.csv',ledger);dump(out/'sex_age_counts.csv',age_rows)
    dump(out/'country_summary.csv',summary);dump(out/'detail_identifiability.csv',detail)
    top=sorted([r for r in ledger if r['sex']=='all' and r['rank_identified'] is not None and r['rank_identified']<=20],key=lambda r:(r['country'],r['rank_identified'],r['code']))
    dump(out/'top20_identified_groups.csv',top)
    result={'schema_version':VERSION,'year':year,'country_scope':raw_scope,'leaf_groups':len(POLICY),'compatible_leaf_groups':sum(v['status']=='compatible' for v in POLICY.values()),'harmonized_rows':len(ledger),'sex_age_rows':len(age_rows),'top20_rows':len(top),'source_audit':audit,'population_age_sum_discrepancies':[{'country':COUNTRIES[c][0],'sex':s,'age_sum_minus_reported_total':v['age_sum_minus_reported_total']} for (c,s),v in pop.items()],'reconciliation':reconciliation,'code_relations':relation_audit(),'new_model_fit':False,'individual_predictions_changed':False,'rates_are_not_individual_probabilities':True,'global_latest_coverage':False,'whole_population_underlying_cause_data':True,'operational_rule':'registered underlying cause; remote infectious aetiology does not reclassify cancer; mixed/ill-defined preserved; residuals use parent complements','no_full_raw_redistribution':True}
    for name,obj in [('results.json',result),('policy.json',{'leaves':POLICY,'parents':PARENTS,'version':VERSION}),('country_coverage.json',raw_scope),('endpoint_relations.json',relation_audit())]:
        (out/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    (out/'SOURCE_MANIFEST.json').write_bytes((source/'manifest.json').read_bytes())
    augment_results(out)
    return result

def augment_results(final: Path):
    """Derived age examples, algebraic identification bounds and model-head links."""
    ledger=list(csv.DictReader((final/'harmonized_ledger.csv').open()))
    ages=list(csv.DictReader((final/'sex_age_counts.csv').open()))
    index={(r['country'],r['sex'],r['code']):r for r in ledger}
    res=[]
    for c in ['RUS','DEU','USA']:
     for component,exactcode,container in [('Alzheimer_G30','1060','1060'),('dementia_F01_F03',None,'1057'),('CKD_N18',None,'1086'),('infective_endocarditis_I33',None,'1068'),('liver_abscess_K750',None,'1080')]:
      current=next(r for r in csv.DictReader((final/'detail_identifiability.csv').open()) if r['country']==c and r['component']==component)
      exact=int(current['deaths']) if current['deaths'] else None
      upper=exact if exact is not None else int(index[c,'all',container]['deaths'])
      lower=exact if exact is not None else 0
      res.append({'country':c,'year':2019,'component':component,'observed_deaths':exact,'identification_lower':lower,'identification_upper':upper,'container_code':container,'interval_type':'algebraic_identification_bound_not_CI','can_use_as_individual_model_coefficient':False})
    dump(final/'identification_bounds.csv',res)
    rate=[]
    for c in ['RUS','DEU']:
     for sex in ['1','2']:
      for cause in ['1067','1069','1052']:
       x=[r for r in ages if r['country']==c and r['sex']==sex and r['code']==cause and r['age'] in ['50-54','55-59']]
       deaths=sum(int(r['deaths']) for r in x);pop=sum(int(r['population']) for r in x)
       rate.append({'country':c,'year':2019,'sex':sex,'age':'50-59','code':cause,'label':POLICY[cause]['label_ru'],'deaths':deaths,'population':pop,'rate_per_100000':deaths/pop*1e5})
    dump(final/'age50_59_examples.csv',rate)
    
    # Map national anchors onto trained heads without assigning subtype coefficients.
    mapping=[]
    for key in sorted({r['code'] for r in csv.DictReader((final/'top20_identified_groups.csv').open())}):
     if key in [str(i) for i in range(1027,1047)]:
      head='malignant_neoplasms';relation='subtype_of_pooled_cancer_head';module='site_specific_oncology';nextstep='individual_records_with_tumour_site_and_death_outcome'
     elif key in ['1065','1067','1068']:
      head='heart_diseases';relation='proper_subset_of_heart_head';module='cardiovascular';nextstep='separate_ICD10_heart_outcomes_and_refit_cause_models'
     elif key=='1066':
      head='heart_diseases_and_other_or_unknown';relation='partial_overlap_I10_I13';module='blood_pressure_cardiorenal';nextstep='align_I10_I11_I12_I13_before_training'
     else:
      head='other_or_unknown';relation='within_heterogeneous_residual';module={'1069':'cerebrovascular','1052':'glycemic_renal','1076':'respiratory_spirometry','1060':'neurodegeneration','1057':'mental_neurodegeneration','1079':'digestive'}.get(key,'disease_specific_review');nextstep='separate_cause_endpoint_and_validate_marker_module'
     mapping.append({'code':key,'cause_ru':POLICY[key]['label_ru'],'current_head':head,'relation':relation,'research_module':module,'published_subtype_coefficients':False,'national_counts_can_set_personal_baseline':False,'next_step':nextstep})
    dump(final/'model_coverage_contract.csv',mapping)
    

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();print(json.dumps(build(a.source,a.out),ensure_ascii=False,indent=2))
