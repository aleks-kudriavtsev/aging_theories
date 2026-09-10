/** Static, browser-only distribution of the exact frozen Python source. */
import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {fileURLToPath} from 'node:url';
const root=path.dirname(fileURLToPath(import.meta.url));
const manifest=JSON.parse(await fs.readFile(path.join(root,'source-files.json'),'utf8'));
const dist=path.join(root,'dist');await fs.rm(dist,{recursive:true,force:true});await fs.mkdir(dist,{recursive:true});
const sources={};
for(const [name,expected] of Object.entries(manifest.files)){
 if(name.includes('..')||name.startsWith('/'))throw Error('Unsafe source path');
 let data;
 if(process.env.LOCAL_SOURCE){data=await fs.readFile(path.join(process.env.LOCAL_SOURCE,name));}
 else {const r=await fetch(`https://raw.githubusercontent.com/${manifest.repository}/${manifest.commit}/${name}`,{signal:AbortSignal.timeout(30000)});if(!r.ok)throw Error('Source retrieval failed: '+name+' '+r.status);data=Buffer.from(await r.arrayBuffer());}
 if(createHash('sha256').update(data).digest('hex')!==expected)throw Error('Source SHA mismatch: '+name);
 sources[name]=data.toString('utf8');
}
const dir='research/mortality/studio16/static/';
let html=sources[dir+'index.html'];
html=html.replace('<meta name="workbench-token" content="__TOKEN__">','<meta name="robots" content="noindex,nofollow"><meta name="description" content="Исследовательское приложение для воспроизведения моделей смертности. Расчёт в браузере без отправки измерений на сервер.">');
html=html.replace('<script src="/app.js" defer>','<script type="module" src="/app.js">');
html=html.replace('0.16</span>','0.16 · Web</span>');
html=html.replace('Сервер работает локально; данные не сохраняются автоматически.','Расчёт выполняется в браузере. Введённые анализы и клинические сведения не отправляются на сервер. Хостинг обрабатывает технические данные загрузки страницы.');
html=html.replace('Проверка файлов…','Загрузка Python в браузер…');
html=html.replace('Не размещайте их в публичном репозитории.','Не размещайте их в публичном репозитории. При закрытии страницы данные не сохраняются приложением.');
html=html.replace('</head>','<link rel="icon" href="/favicon.svg" type="image/svg+xml"></head>');
html=html.replace('</footer>',' <a href="/privacy.html">Обработка данных и область применения</a></footer>');
html=html.replace('</body>','<noscript>Для расчёта необходимо включить JavaScript и поддержку WebAssembly. Данные обрабатываются только в браузере.</noscript></body>');
let js="import {browserFetch} from './browser-api.mjs';\n"+sources[dir+'app.js'].replaceAll('fetch(', 'browserFetch(');
js=js.replace("'X-Workbench-Token':document.querySelector('meta[name=workbench-token]').content", "'X-Client-Mode':'in-browser'");
js=js.replace("body:JSON.stringify(req)","body:selectedMode()==='json'?(comparison?'{\"input\":'+$('jsonInput').value+',\"models\":[\"clinical3\",\"compact4\",\"six6\",\"full8\"]}':$('jsonInput').value):JSON.stringify(req)");
js=js.replace("'Готово. Файл не сохранялся сервером.'", "'Готово. Расчёт выполнен в браузере; данные не отправлялись на сервер.'");
js=js.replace("function collect(){if(selectedMode()==='json')return JSON.parse($('jsonInput').value);", "function collect(){if(selectedMode()==='json')return null;");
const sourceBundle={...manifest,files:{}};
for(const [name,data] of Object.entries(sources))if(name.endsWith('.py')||name.endsWith('.json'))sourceBundle.files[name]={sha256:manifest.files[name],data};
await fs.writeFile(path.join(dist,'python-app.json'),JSON.stringify(sourceBundle));
await fs.writeFile(path.join(dist,'index.html'),html);await fs.writeFile(path.join(dist,'app.js'),js);
await fs.writeFile(path.join(dist,'style.css'),sources[dir+'style.css']);
await fs.writeFile(path.join(dist,'LICENSE.txt'),sources.LICENSE);
for(const name of ['engine-worker.js','browser-api.mjs','privacy.html','favicon.svg'])await fs.copyFile(path.join(root,name),path.join(dist,name));
await fs.writeFile(path.join(dist,'robots.txt'),'User-agent: *\nDisallow: /\n');
if(!process.env.SKIP_VENDOR){
 const v=path.join(dist,'vendor','pyodide-0.28.3');await fs.mkdir(v,{recursive:true});
 for(const name of ['pyodide.js','pyodide.mjs','pyodide.asm.js','pyodide.asm.wasm','python_stdlib.zip','pyodide-lock.json'])await fs.copyFile(path.join(root,'node_modules','pyodide',name),path.join(v,name));
 for(const name of ['LICENSE','LICENSE.txt','LICENSE.md'])try{await fs.copyFile(path.join(root,'node_modules','pyodide',name),path.join(v,name));}catch(e){if(e.code!=='ENOENT')throw e;}
}
await fs.writeFile(path.join(dist,'deployment-info.json'),JSON.stringify({application:'Mortality Studio 0.16 Web',source_commit:manifest.commit,model_files:Object.fromEntries(Object.entries(manifest.files).filter(([n])=>n.endsWith('.json'))),runtime:'Pyodide 0.28.3 / CPython in Web Worker',input_processing:'browser_only',backend_api:false,analytics:false,source_files:Object.keys(sourceBundle.files).length},null,2));
console.log('STATIC_BUILD_OK: verified sources, unchanged models, browser-only processing.');
