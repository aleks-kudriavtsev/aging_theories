/* No user code or request-body logging. Only fixed Python functions. */
'use strict';
let dispatch=null;
const ready=(async()=>{
 importScripts('/vendor/pyodide-0.28.3/pyodide.js');
 const py=await loadPyodide({indexURL:'/vendor/pyodide-0.28.3/',stdout:()=>{},stderr:()=>{},fullStdLib:false});
 const response=await fetch('/python-app.json',{cache:'no-cache'});
 if(!response.ok)throw Error('application_bundle_unavailable');
 const bundle=await response.json();
 for(const [name,row] of Object.entries(bundle.files)){
  if(name.includes('..')||name.startsWith('/')||!name.startsWith('research/'))throw Error('invalid_bundle_path');
  const bytes=new TextEncoder().encode(row.data);
  const hash=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',bytes)),x=>x.toString(16).padStart(2,'0')).join('');
  if(hash!==row.sha256)throw Error('source_integrity_failure');
  const p='/application/'+name;py.FS.mkdirTree(p.slice(0,p.lastIndexOf('/')));py.FS.writeFile(p,bytes);
 }
 py.runPython(`
import sys, json
sys.path.insert(0, '/application')
from urllib.parse import urlsplit, parse_qs
from research.mortality.studio16.registry import catalog
from research.mortality.studio16.engine import run, compare, example
from research.mortality.studio16.report import render
from research.mortality.workbench.runtime import strict_json

def _browser_dispatch(path, method, raw):
    try:
        parsed = urlsplit(path)
        if method == 'GET' and parsed.path == '/api/catalog':
            data = catalog()
        elif method == 'GET' and path == '/api/health':
            catalog()
            data = {'status':'ready_research_only','clinical_use_ready':False,'processing':'browser_only'}
        elif method == 'GET' and parsed.path == '/api/example':
            args = parse_qs(parsed.query, strict_parsing=True)
            if set(args) != {'model_id'} or len(args['model_id']) != 1: raise ValueError('model_required')
            data = example(args['model_id'][0])
        elif method == 'POST' and path in ('/api/predict', '/api/compare'):
            p = strict_json(raw)
            if path == '/api/compare':
                if not isinstance(p, dict) or set(p) != {'input','models'}: raise ValueError('invalid_envelope')
                result = compare(p['input'], p['models'])
            else:
                result = run(p)
            result['execution'] = {'location':'browser','measurements_sent_to_server':False,'source_commit':'5f51320b7556b8afce2c311797ed6b4cbac4d0a6'}
            data = {'result':result,'report_html':render(result)}
            return json.dumps({'status':422 if result['status']=='blocked' else 200,'data':data}, ensure_ascii=False,allow_nan=False)
        else:
            raise ValueError('unsupported_local_command')
        return json.dumps({'status':200,'data':data},ensure_ascii=False,allow_nan=False)
    except (ValueError,KeyError,TypeError,ArithmeticError,RecursionError,OSError):
        return json.dumps({'status':400,'data':{'error':'invalid_input_or_integrity_failure','clinical_use_ready':False}})
`);
 dispatch=py.globals.get('_browser_dispatch');
 globalThis.fetch=()=>Promise.reject(new Error('Network disabled after initialization'));
 globalThis.XMLHttpRequest=undefined;globalThis.WebSocket=undefined;
 return true;
})();
ready.catch(e=>self.postMessage({kind:'startup_error',message:String(e.message)}));
self.onmessage=async event=>{
 const {id,path,method,body}=event.data||{};
 try{
  await ready;
  if(!Number.isSafeInteger(id)||!['GET','POST'].includes(method)||typeof path!=='string'||typeof body!=='string'||body.length>1048576)throw Error('invalid_command');
  self.postMessage({id,payload:JSON.parse(dispatch(path,method,body))});
 }catch(e){self.postMessage({id,payload:{status:503,data:{error:'Браузерный расчёт недоступен. Обновите страницу в современном браузере; данные не отправлены.',clinical_use_ready:false}}});}
};
