"""Self-contained HTML export, all variable text escaped. No external resources."""
from html import escape
import json
from pathlib import Path
from .registry import LABELS


def render(result):
    def text(x):return escape(str(x))
    sections=[]
    models=result.get('results',[result])
    for model in models:
        parts=['<section><h2>'+text(model.get('model_title','Расчёт заблокирован'))+'</h2>']
        for row in model.get('probabilities') or []:
            parts+=['<h3>Горизонт: '+text(row['horizon_years'])+' лет</h3><table><tr><th>Исход</th><th>Вероятность</th></tr>']
            for k,v in [('Выживание',row['survival']),('Все причины смерти',row['all_cause_death_probability'])]+[(LABELS.get(k,k),v) for k,v in row['cause_specific_cif'].items()]:
                parts+=['<tr><td>'+text(k)+'</td><td>'+text(f'{v:.3%}')+'</td></tr>']
            parts+=['</table><p>Все причины — сумма причин смерти, не дополнительная категория.</p>']
        parts+=['<h3>Проверка на исследовательских данных</h3>']
        for v in model.get('validation',[]):
            parts+=['<p>'+text(v.get('cohort',''))+'; горизонт '+text(v['horizon_years'])+'; '+text(v.get('caveat',''))+'</p>']
            if 'all_cause_auc' in v:parts+=['<p>AUC на проверочной когорте: '+text(f"{v['all_cause_auc']:.4f}")+'. Это не точность текущего индивидуального прогноза.</p>']
        if model.get('clipped_features'):parts+=['<p><strong>Ограничены обучающим диапазоном: '+text(', '.join(model['clipped_features']))+'</strong></p>']
        parts+=['<p>SHA модели: '+text(model.get('model_sha256','—'))+'</p></section>']
        sections+=parts
    return '<!doctype html><html lang="ru"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Mortality Studio 0.16 — исследовательский отчёт</title><style>body{font:16px/1.5 system-ui;max-width:1000px;margin:40px auto;padding:0 24px;color:#18232e}table{border-collapse:collapse;width:100%}th,td{padding:8px;text-align:left;border-bottom:1px solid #ddd}section{break-inside:avoid;margin:24px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}aside{border:2px solid #9c6425;padding:14px}h1{font-size:28px}</style><h1>Mortality Studio 0.16</h1><aside>'+text(result.get('disclaimer','Только исследовательское использование.'))+'</aside>'+''.join(sections)+'<details><summary>Полные результаты и происхождение</summary><pre>'+text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))+'</pre></details><p>Файл содержит введённые измерения и результаты. Храните его как исследовательские данные; не публикуйте автоматически.</p></html>'


def save(result,out):
    out=Path(out)
    if out.exists() or out.is_symlink():raise ValueError('Choose a new report directory')
    # Serialize before creating output; failures cannot leave a misleading report.
    raw=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)
    html=render(result);out.mkdir(parents=True)
    (out/'result.json').write_text(raw+'\n',encoding='utf-8')
    (out/'report.html').write_text(html,encoding='utf-8')
