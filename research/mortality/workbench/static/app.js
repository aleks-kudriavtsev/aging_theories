'use strict';
const $ = id => document.getElementById(id);
let spec, latestResult = null, latestHTML = null, inputRevision = 0;
const clinicalNames = {smoking:'Курение',bp_treatment:'Антигипертензивное лечение',diabetes_history:'Диабет в анамнезе',cvd_history:'ССЗ в анамнезе',cancer_history:'Рак в анамнезе'};
const demoValues = {bmi:26,sbp:125,total_cholesterol:195,hdl:55,hba1c:5.4,creatinine:0.9,uacr:8,albumin:43,rdw:12.8,wbc:6.5,crp:1.5,ntprobnp:60,hs_ctnt:6,cystatin_c:0.8};
function el(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n;}
function option(select,value,text){const o=el('option',text);o.value=value;select.append(o);}
function labelled(text,node){const l=el('label',text);l.append(node);return l;}
function invalidateResult(){
  inputRevision++;
  if(!latestResult && $('output').hidden)return;
  latestResult=null;latestHTML=null;
  $('risk-cards').replaceChildren();$('protocol').textContent='';
  $('result-status').textContent='Требуется новый расчёт';
  $('output-context').textContent='Ввод изменён. Предыдущий результат больше не относится к текущим данным.';
  $('output-message').replaceChildren();
  for(const id of ['download-html','download-json'])$(id).disabled=true;
}
function buildFields(){
  const saved={};
  for(const node of document.querySelectorAll('#clinical-fields select, #measurements input, #measurements select'))saved[node.id]=node.value;
  const model=$('model').value; $('clinical').hidden=model==='M0_age_sex'; $('clinical-fields').replaceChildren();
  for(const key of spec.clinical_fields){const s=el('select');s.id='c-'+key;option(s,'','Выберите');
    if(key==='smoking'){for(const [v,t] of [['never','Никогда'],['former','В прошлом'],['current','Сейчас']])option(s,v,t);}
    else {option(s,'false','Нет');option(s,'true','Да');} $('clinical-fields').append(labelled(clinicalNames[key],s));}
  $('measurements').replaceChildren();
  for(const key of spec.required_measurements[model]){
    const f=spec.measurements[key],row=el('div',undefined,'measurement');row.append(el('div',f.label,'caption'));
    const v=el('input');v.type='number';v.step='any';v.id='v-'+key;v.setAttribute('aria-label',f.label);row.append(v);
    const u=el('select');u.id='u-'+key;u.setAttribute('aria-label','Единицы: '+f.label);for(const unit of Object.keys(f.conversions))option(u,unit,unit);row.append(u);
    const meta=el('div',undefined,'meta'),m=el('select');m.id='m-'+key;option(m,f.matrix,f.matrix+' — требование модели');option(m,'other','Другой / неизвестен');
    meta.append(labelled('Материал / источник измерения',m));
    const a=el('select');a.id='a-'+key;option(a,'not_reported','Не подтверждён');if(f.required_assay)option(a,f.required_assay,f.required_assay);else option(a,'method_recorded_separately','Метод описан в исходном протоколе');
    meta.append(labelled(f.required_assay?'Подтвердите метод / стандартизацию':'Метод (переносимость не установлена)',a));row.append(meta);$('measurements').append(row);
  }
  for(const [id,value] of Object.entries(saved)){const node=$(id);if(node)node.value=value;}
}
function payload(){
  const model=$('model').value,clinical={},measurements={};
  if(model!=='M0_age_sex')for(const k of spec.clinical_fields){const v=$('c-'+k).value;clinical[k]=k==='smoking'?(v||null):v===''?null:v==='true';}
  for(const key of spec.required_measurements[model]){const text=$('v-'+key).value;measurements[key]={value:text===''?null:Number(text),unit:$('u-'+key).value,matrix:$('m-'+key).value,assay:$('a-'+key).value};}
  return {schema_version:'0.6.0',research_only:$('ack').checked,acknowledge_limitations:$('ack').checked,
    dataset_context:$('context').value,country:$('country').value,endpoint:'all_cause_death',model,
    horizons_years:[1,5,10],age_years:$('age').value===''?null:Number($('age').value),sex:$('sex').value||null,clinical,measurements};
}
function showError(text){$('output').hidden=false;$('output-message').replaceChildren(el('p',text,'error'));$('risk-cards').replaceChildren();$('result-status').textContent='Ошибка';$('output-context').textContent='';$('protocol').textContent='';latestHTML=null;latestResult=null;for(const id of ['download-html','download-json'])$(id).disabled=true;}
async function submit(raw){
  const revision=inputRevision;
  latestHTML=null;latestResult=null;
  for(const id of ['download-html','download-json'])$(id).disabled=true;
  for(const id of ['calculate','calculate-json'])$(id).disabled=true;
  try{const response=await fetch('/api/predict',{method:'POST',headers:{'Content-Type':'application/json','X-Workbench-Token':document.querySelector('meta[name=workbench-token]').content},body:raw});
    const data=await response.json();
    if(revision!==inputRevision){showError('Данные изменены во время запроса. Выполните новый расчёт.');return;}
    if(!data.result)throw Error(data.error||'Ошибка ответа');
    const r=data.result;latestResult=r;latestHTML=data.report_html;$('output').hidden=false;$('result-status').textContent=r.status==='blocked'?'Расчёт заблокирован':'RESEARCH ONLY';
    $('output-context').textContent=r.declaration==='synthetic_not_a_person'?'Синтетический пример: результат не относится к человеку.':r.declaration==='user_declared_historical_research'?'Заявленные исторические исследовательские данные США.':'Область применения не подтверждена.';
    $('risk-cards').replaceChildren();$('output-message').replaceChildren();
    if(r.probabilities){for(const x of r.probabilities){const c=el('div',undefined,'risk');c.append(el('div',x.horizon_years+' лет'),el('strong',(x.probability*100).toFixed(2)+'%'),el('small','Выход исторической модели, не клинический прогноз'));$('risk-cards').append(c);}
      $('output-message').append(el('p','Нет установленного индивидуального доверительного интервала и распределения по причинам смерти. Никаких рекомендаций по лечению не формируется.','muted'));
    }else $('output-message').append(el('p','Вероятности не рассчитаны. Исправьте поля, перечисленные ниже.','error'));
    $('protocol').textContent=JSON.stringify(r,null,2);for(const id of ['download-html','download-json'])$(id).disabled=false;$('output').scrollIntoView({behavior:'smooth'});
  }catch(e){showError('Не удалось выполнить расчёт: '+e.message);}finally{for(const id of ['calculate','calculate-json'])$(id).disabled=false;}
}
function download(name,contents,type){const url=URL.createObjectURL(new Blob([contents],{type}));const a=el('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
async function init(){
 spec=await (await fetch('/api/schema')).json();for(const [k,v] of Object.entries(spec.models))option($('model'),k,v);$('model').value='M1_routine';buildFields();
 $('model').addEventListener('change',buildFields);
 $('profile').addEventListener('input',invalidateResult);
 $('profile').addEventListener('change',invalidateResult);
 $('json-input').addEventListener('input',invalidateResult);
 $('demo').addEventListener('click',()=>{invalidateResult();$('age').value=55;$('sex').value='female';$('country').value='US';$('context').value='synthetic_example';$('ack').checked=true;
  for(const k of spec.clinical_fields)$('c-'+k).value=k==='smoking'?'never':'false';
  for(const k of spec.required_measurements[$('model').value]){$('v-'+k).value=demoValues[k];$('u-'+k).value=spec.measurements[k].unit;$('m-'+k).value=spec.measurements[k].matrix;$('a-'+k).value=spec.measurements[k].required_assay||'not_reported';}
  $('json-input').value=JSON.stringify(payload(),null,2);
 });
 $('profile').addEventListener('submit',e=>{e.preventDefault();submit(JSON.stringify(payload()));});
 $('to-json').addEventListener('click',()=>{invalidateResult();$('json-input').value=JSON.stringify(payload(),null,2);});
 $('calculate-json').addEventListener('click',()=>submit($('json-input').value));
 $('upload').addEventListener('change',async()=>{invalidateResult();const f=$('upload').files[0];if(!f)return;if(f.size>1048576)return showError('Допускается JSON не более 1 MiB.');$('json-input').value=await f.text();});
 $('download-html').addEventListener('click',()=>{if(latestHTML)download('research_report.html',latestHTML,'text/html');});
 $('download-json').addEventListener('click',()=>{if(latestResult)download('research_report.json',JSON.stringify(latestResult,null,2),'application/json');});
}
init().catch(e=>showError('Не удалось загрузить интерфейс: '+e.message));
