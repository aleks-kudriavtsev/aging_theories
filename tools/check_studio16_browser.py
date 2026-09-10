import sys,threading,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import argparse,os
parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
OUT=args.out
from research.mortality.studio16.server import create_server
from playwright.sync_api import sync_playwright, expect
s=create_server(0);t=threading.Thread(target=s.serve_forever,daemon=True);t.start()
url=f'http://127.0.0.1:{s.server_port}'
checks=[]
try:
 with sync_playwright() as p:
  browser=p.chromium.launch(headless=True, **({'executable_path':os.environ['CHROMIUM_EXECUTABLE']} if os.environ.get('CHROMIUM_EXECUTABLE') else {}))
  page=browser.new_page(viewport={'width':1460,'height':1080})
  errors=[];page.on('pageerror',lambda err:errors.append(str(err)))
  page.goto(url,wait_until='networkidle',timeout=20000)
  expect(page.locator('#health')).to_contain_text('проверены')
  checks.append('real_local_navigation')
  page.select_option('#model','compact4');page.click('#example')
  expect(page.locator('#age')).to_have_value('55')
  page.click('#predict');expect(page.locator('#downloadJSON')).to_be_enabled()
  assert page.locator('.risk').count()==2
  checks.append('form_to_actual_API_two_horizons')
  page.screenshot(path=str(OUT/'studio16_desktop.png'),full_page=True)
  page.fill('#age','60');assert page.locator('#downloadJSON').is_disabled();assert page.locator('.risk').count()==0
  checks.append('result_invalidated_on_edit')
  page.select_option('#country','RU');page.click('#predict');expect(page.locator('#message')).to_contain_text('не выполнены')
  assert page.locator('.risk').count()==0
  checks.append('Russia_blocked_no_probabilities')
  page.select_option('#model','full8');page.click('#example');expect(page.locator('#age')).to_have_value('55')
  page.click('#compare');expect(page.locator('.prediction')).to_have_count(4)
  assert page.locator('.risk').count()==8
  checks.append('four_panel_comparison_actual_API')
  with page.expect_download() as download:
   page.click('#downloadHTML')
  download.value.save_as(str(OUT/'browser_download.html'))
  checks.append('download_HTML')
  page.check('input[name=mode][value=json]')
  assert page.locator('#profile').is_hidden()
  payload=json.loads(page.input_value('#jsonInput'));payload['model_id']='compact4';payload['measurements']={k:v for k,v in payload['measurements'].items() if k in ['bmi','sbp','hba1c','creatinine','uacr','albumin']}
  page.fill('#jsonInput',json.dumps(payload));page.click('#predict');expect(page.locator('.prediction')).to_have_count(1)
  assert page.locator('.prediction h3').inner_text().startswith('Компактная')
  assert 'Компактная' in page.locator('#evidence').inner_text()
  checks.append('JSON_model_identity_and_evidence_actual_API')
  page.set_viewport_size({'width':390,'height':844});page.screenshot(path=str(OUT/'studio16_mobile.png'),full_page=True)
  assert page.evaluate('() => document.documentElement.scrollWidth <= innerWidth')
  checks.append('mobile_no_horizontal_overflow')
  assert not errors,errors
  checks.append('no_javascript_errors')
  browser.close()
except Exception as e:
 (OUT/'browser_status.json').write_text(json.dumps({'passed':checks,'error':str(e),'real_browser_end_to_end_complete':False},indent=2));raise
finally:s.shutdown();s.server_close();t.join()
(OUT/'browser_status.json').write_text(json.dumps({'passed':checks,'count':len(checks),'real_browser_end_to_end_complete':True,'mocked_network':False},indent=2))
print(checks)
