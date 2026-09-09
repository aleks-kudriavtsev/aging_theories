"""Escaped, standalone HTML report and atomic local output. No external assets."""
from __future__ import annotations
import html
import json
from pathlib import Path
import shutil
import tempfile
from .runtime import DISCLAIMER

STYLE = '''body{font:16px/1.55 system-ui,sans-serif;background:#f3f5f7;color:#182631;margin:0}
main{max-width:1120px;margin:40px auto;padding:36px;background:white;border-radius:18px}
h1{font-size:30px;line-height:1.2;margin:8px 0 20px}h2{margin-top:32px;font-size:21px}
small,.muted{color:#536676}.notice{padding:16px 20px;background:#fff5dd;border-left:4px solid #aa7712;border-radius:6px}
.grid{display:flex;gap:18px;flex-wrap:wrap;margin:24px 0}.card{min-width:180px;flex:1;padding:18px;background:#edf3f5;border-radius:12px}
.value{font-size:30px;font-weight:650}table{width:100%;border-collapse:collapse;font-size:14px}th,td{text-align:left;padding:9px;border-bottom:1px solid #dce3e8;vertical-align:top}
pre,code{font-size:12px;white-space:pre-wrap;overflow-wrap:anywhere}.scroll{overflow:auto}.blocked{font-weight:bold}details{margin:14px 0}
@media(max-width:700px){main{margin:0;padding:18px;border-radius:0}.card{min-width:130px}}@media print{body{background:white}main{margin:0;padding:0}.card{break-inside:avoid}details{display:block}}'''


def esc(value):return html.escape(str(value),quote=True)


def render(result):
    status=esc(result.get('status','invalid'))
    if result.get('probabilities') is not None:
        cards=''.join('<div class="card"><div>'+esc(r['horizon_years'])+' лет</div><div class="value">'+
            f"{r['probability']*100:.2f}%"+'</div><small>Историческая модель, не клинический прогноз</small></div>' for r in result['probabilities'])
        body='<div class="grid">'+cards+'</div>'
        rows=''.join('<tr>'+''.join('<td>'+esc(r.get(k,''))+'</td>' for k in ['field','input_value','input_unit','canonical_value','canonical_unit'])+'</tr>' for r in result.get('measurement_trace',[]))
        body+='<h2>Проверка единиц</h2><div class="scroll"><table><tr><th>Показатель</th><th>Вход</th><th>Единица входа</th><th>После преобразования</th><th>Единица модели</th></tr>'+rows+'</table></div>'
        body+='<p>Импутация отсутствующих данных: <strong>не применялась</strong>. Разделение смертности по причинам: <strong>не рассчитано</strong>. Индивидуальный доверительный интервал: <strong>не установлен</strong>.</p>'
    else:
        body='<p class="blocked">Расчёт заблокирован. Вероятности не сформированы.</p>'
    body+='<h2>Ошибки и ограничения</h2><pre>'+esc(json.dumps({'errors':result.get('errors',[]),'warnings':result.get('warnings',[])},ensure_ascii=False,indent=2))+'</pre>'
    body+='<details><summary>Полный протокол вычисления и преобразований</summary><pre>'+esc(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))+'</pre></details>'
    return ('<!doctype html><html lang="ru"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; style-src \'unsafe-inline\'; base-uri \'none\'; form-action \'none\'">'
            '<title>Исследовательский отчёт — Mortality Workbench</title><style>'+STYLE+'</style></head><body><main>'
            '<small>MORTALITY WORKBENCH · 0.7.0 · RESEARCH ONLY</small><h1>Протокол исследовательского расчёта</h1>'
            '<div class="notice">'+esc(DISCLAIMER)+'</div>'
            '<p><b>Статус:</b> '+status+'<br><b>Модель:</b> '+esc(result.get('model','—'))+
            '<br><b>Контекст:</b> '+esc(result.get('declaration','Не установлен'))+'</p>'+body+
            '<h2>Происхождение модели</h2><code>База кода: '+esc(result.get('model_code_base_commit',''))+'<br>SHA-256: '+esc(result.get('model_sha256',''))+'</code>'
            '<p class="muted">Нет биологического возраста, рекомендаций по лечению или автоматического распределения риска по долям ВОЗ. Сохранённый отчёт содержит введённые измерения; храните его локально.</p>'
            '</main></body></html>')


def write_package(result, out):
    """Write to a NEW directory; never overwrite prior data or follow symlinks."""
    out=Path(out)
    if out.exists() or out.is_symlink():raise ValueError('Каталог результата уже существует: выберите новый')
    out.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(tempfile.mkdtemp(prefix='.mortality-',dir=out.parent))
    try:
        (tmp/'report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        (tmp/'report.html').write_text(render(result),encoding='utf-8')
        # Reserve destination exclusively. Copying into this new directory does
        # not modify another report; remove it on any exceptional write failure.
        out.mkdir()
        try:
            for p in tmp.iterdir():shutil.move(str(p),str(out/p.name))
        except Exception:
            shutil.rmtree(out);raise
    finally:shutil.rmtree(tmp,ignore_errors=True)
    return out
