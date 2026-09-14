"""Loopback-only research UI. Never use this development server publicly."""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import json
import secrets
from urllib.parse import urlsplit, parse_qs
from .registry import catalog
from .engine import run, compare, example
from .report import render
from ..workbench.runtime import MAX_BODY, strict_json
STATIC=Path(__file__).with_name('static')


def doctor():
    try:
        c=catalog()
        for name in ('index.html','app.js','style.css'):
            if not (STATIC/name).is_file() or not (STATIC/name).stat().st_size:raise ValueError('Missing asset')
    except (OSError,ValueError,KeyError,TypeError):
        return {'version':'0.16.0','status':'not_ready','clinical_use_ready':False}
    return {'version':'0.16.0','status':'ready_research_only','registered_views':len(c['models']),
        'frozen_artifacts_verified':3,'clinical_use_ready':False,'Russia_calibrated':False,
        'Germany_calibrated':False,'new_model_fit':False,'network_binding':'127.0.0.1',
        'request_storage':False,'external_requests_by_server':False,'legacy_cli_preserved':True}


def create_server(port=8766):
    if type(port) is not int or not 0<=port<=65535:raise ValueError('Invalid port')
    if doctor()['status']!='ready_research_only':raise ValueError('Incomplete release')
    token=secrets.token_urlsafe(32)
    class Handler(BaseHTTPRequestHandler):
        server_version='MortalityStudio/0.16';sys_version=''
        def setup(self):super().setup();self.connection.settimeout(5)
        def log_message(self,*args):pass
        def host_ok(self):return self.headers.get_all('Host')==[f'127.0.0.1:{self.server.server_port}']
        def reply(self,code,data,ctype='application/json; charset=utf-8'):
            raw=data.encode('utf-8') if isinstance(data,str) else data
            self.send_response(code);self.send_header('Content-Type',ctype);self.send_header('Content-Length',str(len(raw)))
            for key,value in [('Cache-Control','no-store'),('X-Content-Type-Options','nosniff'),('Referrer-Policy','no-referrer'),('X-Frame-Options','DENY'),('Content-Security-Policy',"default-src 'self'; connect-src 'self'; script-src 'self'; style-src 'self'; img-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'")]:self.send_header(key,value)
            self.end_headers();self.wfile.write(raw)
        def data(self,code,data):self.reply(code,json.dumps(data,ensure_ascii=False,allow_nan=False))
        def error(self,code,message):self.data(code,{'error':message,'clinical_use_ready':False})
        def do_GET(self):
            if not self.host_ok():return self.error(403,'host_not_allowed')
            url=urlsplit(self.path)
            if url.path=='/api/health':
                d=doctor();return self.data(200 if d['status']=='ready_research_only' else 503,d)
            try:
                if url.path=='/api/catalog':return self.data(200,catalog())
                if url.path=='/api/example':
                    args=parse_qs(url.query,strict_parsing=True)
                    if set(args)!={'model_id'} or len(args['model_id'])!=1:raise ValueError('Model required')
                    return self.data(200,example(args['model_id'][0]))
                routes={'/':('index.html','text/html; charset=utf-8'),'/app.js':('app.js','application/javascript; charset=utf-8'),'/style.css':('style.css','text/css; charset=utf-8')}
                if self.path not in routes:return self.error(404,'not_found')
                name,ctype=routes[self.path];content=(STATIC/name).read_text(encoding='utf-8')
                self.reply(200,content.replace('__TOKEN__',token) if name=='index.html' else content,ctype)
            except (OSError,ValueError,KeyError,TypeError):return self.error(400,'invalid_request_or_artifact')
        def do_POST(self):
            if not self.host_ok():return self.error(403,'host_not_allowed')
            if self.path not in ('/api/predict','/api/compare'):return self.error(404,'not_found')
            origins=self.headers.get_all('Origin')
            if origins is not None and origins!=[f'http://127.0.0.1:{self.server.server_port}']:return self.error(403,'origin_not_allowed')
            tokens=self.headers.get_all('X-Workbench-Token') or []
            if len(tokens)!=1 or not secrets.compare_digest(tokens[0].encode('utf-8'),token.encode('ascii')):return self.error(403,'request_token_required')
            if self.headers.get_all('Transfer-Encoding'):return self.error(400,'transfer_encoding_not_supported')
            lengths=self.headers.get_all('Content-Length') or []
            if len(lengths)!=1:return self.error(411,'one_content_length_required')
            try:n=int(lengths[0])
            except ValueError:return self.error(400,'invalid_length')
            if n<=0 or n>MAX_BODY:return self.error(413,'payload_too_large_or_empty')
            if (self.headers.get('Content-Type','').split(';')[0].strip()!='application/json'
                or len(self.headers.get_all('Content-Type') or [])!=1):return self.error(415,'application_json_required')
            if doctor()['status']!='ready_research_only':return self.error(503,'not_ready')
            try:
                raw=self.rfile.read(n)
                if len(raw)!=n:raise ValueError('Incomplete body')
                p=strict_json(raw)
                if self.path=='/api/compare':
                    if not isinstance(p,dict) or set(p)!={'input','models'}:raise ValueError('Invalid comparison envelope')
                    result=compare(p['input'],p['models'])
                else:result=run(p)
                return self.data(422 if result['status']=='blocked' else 200,{'result':result,'report_html':render(result)})
            except (ValueError,KeyError,TypeError,RecursionError,OverflowError,ArithmeticError,TimeoutError,OSError):return self.error(400,'invalid_input_or_integrity_failure')
    server=ThreadingHTTPServer(('127.0.0.1',port),Handler);server.daemon_threads=True
    return server


def serve(port=8766):
    with create_server(port) as server:
        print(f'Исследовательский интерфейс: http://127.0.0.1:{server.server_port}',flush=True)
        print('Данные не сохраняются сервером. Только локальный доступ. Ctrl+C — остановить.',flush=True)
        try:server.serve_forever()
        except KeyboardInterrupt:pass
