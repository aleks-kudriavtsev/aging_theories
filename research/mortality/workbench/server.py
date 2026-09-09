"""Loopback-only development UI, no file server, data persistence or access log.

Not a production server. Host/Origin checks, a per-process POST token and strict
JSON reduce accidental cross-origin use; they do not replace a security review.
"""
from __future__ import annotations
from functools import partial
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import secrets
from .runtime import MAX_BODY, run, schema, strict_json, load_bundle
from .report import render

STATIC=Path(__file__).with_name('static')


def check_assets():
    """Fail early when a copied/installed distribution lacks UI resources."""
    for name in ('index.html', 'app.js', 'style.css'):
        path = STATIC / name
        if not path.is_file() or not path.stat().st_size:
            raise ValueError('Missing required UI asset: ' + name)


def readiness():
    try:
        load_bundle()
        check_assets()
    except (ValueError, OSError, KeyError, TypeError):
        return {'application_version': '0.7.0', 'status': 'not_ready',
                'clinical_use_ready': False}
    return {'application_version': '0.7.0', 'status': 'ready_research_only',
            'clinical_use_ready': False}


def create_server(port=8765):
    if type(port) is not int or not 0<=port<=65535:raise ValueError('Invalid local port')
    check_assets()
    token=secrets.token_urlsafe(32)
    class Handler(BaseHTTPRequestHandler):
        server_version='MortalityWorkbench/0.7'
        sys_version=''
        def setup(self):
            super().setup();self.connection.settimeout(5)
        def log_message(self,*args):pass  # no submitted data in stdout or access logs
        def check_host(self):
            expected=f'127.0.0.1:{self.server.server_port}'
            return self.headers.get_all('Host')==[expected]
        def reply(self,code,data,ctype='application/json; charset=utf-8'):
            raw=data.encode('utf-8') if isinstance(data,str) else data
            self.send_response(code);self.send_header('Content-Type',ctype)
            self.send_header('Content-Length',str(len(raw)))
            self.send_header('Cache-Control','no-store')
            self.send_header('X-Content-Type-Options','nosniff')
            self.send_header('Referrer-Policy','no-referrer')
            self.send_header('X-Frame-Options','DENY')
            self.send_header('Content-Security-Policy',"default-src 'self'; connect-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'")
            self.end_headers();self.wfile.write(raw)
        def error(self,code,message):
            self.reply(code,json.dumps({'error':message},ensure_ascii=False))
        def do_GET(self):
            if not self.check_host():return self.error(403,'host_not_allowed')
            if self.path=='/api/health':
                status=readiness()
                return self.reply(200 if status['status']=='ready_research_only' else 503,json.dumps(status))
            if self.path=='/api/schema':return self.reply(200,json.dumps(schema(),ensure_ascii=False))
            routes={'/':('index.html','text/html; charset=utf-8'),
                    '/app.js':('app.js','application/javascript; charset=utf-8'),
                    '/style.css':('style.css','text/css; charset=utf-8')}
            if self.path not in routes:return self.error(404,'not_found')
            name,ctype=routes[self.path];content=(STATIC/name).read_text(encoding='utf-8')
            if name=='index.html':content=content.replace('__TOKEN__',token)
            self.reply(200,content,ctype)
        def do_POST(self):
            if not self.check_host():return self.error(403,'host_not_allowed')
            if self.path!='/api/predict':return self.error(404,'not_found')
            origin=self.headers.get('Origin')
            if origin is not None and origin!=f'http://127.0.0.1:{self.server.server_port}':
                return self.error(403,'origin_not_allowed')
            if not secrets.compare_digest(self.headers.get('X-Workbench-Token','').encode('utf-8'),token.encode('ascii')):
                return self.error(403,'request_token_required')
            if self.headers.get('Transfer-Encoding'):return self.error(400,'transfer_encoding_not_supported')
            lengths=self.headers.get_all('Content-Length') or []
            if len(lengths)!=1:return self.error(411,'one_content_length_required')
            try:n=int(lengths[0])
            except ValueError:return self.error(400,'invalid_length')
            if n<=0 or n>MAX_BODY:return self.error(413,'payload_too_large_or_empty')
            if self.headers.get('Content-Type','').split(';')[0].strip()!='application/json':
                return self.error(415,'application_json_required')
            if readiness()['status']!='ready_research_only':
                return self.error(503,'application_not_ready')
            try:
                raw=self.rfile.read(n)
                if len(raw)!=n:raise ValueError('incomplete_body')
                payload=strict_json(raw);result=run(payload)
            except (ValueError,KeyError,TypeError,RecursionError,OverflowError,TimeoutError,OSError):
                return self.error(400,'invalid_input_or_model_integrity_failure')
            response={'result':result,'report_html':render(result)}
            self.reply(422 if result['status']=='blocked' else 200,json.dumps(response,ensure_ascii=False,allow_nan=False))
    server=ThreadingHTTPServer(('127.0.0.1',port),Handler)
    server.daemon_threads=True
    return server


def serve(port=8765):
    from .runtime import load_bundle
    load_bundle()
    with create_server(port) as server:
        print(f'Локальный интерфейс: http://127.0.0.1:{server.server_port}',flush=True)
        print('Только исследование. Данные не записываются сервером. Остановка: Ctrl+C.',flush=True)
        try:server.serve_forever()
        except KeyboardInterrupt:pass
