/** API-shaped in-memory messages; does not perform HTTP requests. */
const worker=new Worker('/engine-worker.js');
const pending=new Map();let id=0;
worker.addEventListener('message',event=>{const item=pending.get(event.data.id);if(!item)return;pending.delete(event.data.id);clearTimeout(item.timer);item.cleanup();item.resolve(event.data.payload);});
worker.addEventListener('error',()=>{for(const item of pending.values()){clearTimeout(item.timer);item.cleanup();item.reject(new Error('Не удалось загрузить Python в браузер. Проверьте соединение и обновите страницу.'));}pending.clear();});
export async function browserFetch(path,options={}){
 const method=options.method||'GET';const body=options.body||'';
 if(!['/api/catalog','/api/health','/api/predict','/api/compare'].includes(path)&&!path.startsWith('/api/example?'))throw Error('Unsupported local route');
 const number=++id;
 const payload=await new Promise((resolve,reject)=>{
  const abort=()=>{const item=pending.get(number);if(!item)return;clearTimeout(item.timer);pending.delete(number);item.cleanup();reject(new DOMException('Input changed','AbortError'));};
  if(options.signal?.aborted){reject(new DOMException('Input changed','AbortError'));return;}
  const cleanup=()=>options.signal?.removeEventListener('abort',abort);
  const timer=setTimeout(()=>{pending.delete(number);cleanup();reject(new Error('Не завершилась загрузка вычислительного модуля. Обновите страницу.'));},120000);
  pending.set(number,{resolve,reject,timer,cleanup});options.signal?.addEventListener('abort',abort,{once:true});
  worker.postMessage({id:number,path,method,body});
 });
 return {ok:payload.status>=200&&payload.status<300,status:payload.status,json:async()=>payload.data};
}
