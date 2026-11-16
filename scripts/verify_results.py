from pathlib import Path
from fractions import Fraction as F
import argparse,json,math,subprocess,hashlib
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads((ROOT/p).read_text())
def poly(p):return {tuple(t['powers']):int(t['coefficient']) for t in load(p)}
def sub(a,b):return {k:a.get(k,0)-b.get(k,0) for k in a.keys()|b.keys() if a.get(k,0)!=b.get(k,0)}
def norm(w):
 d={};return tuple(d.setdefault(a,len(d)) for a in w)
def canonical(w):return min(norm(v[i:]+v[:i]) for v in [w,w[::-1]] for i in range(len(w)))
def graph(w):
 p={a:[i for i,x in enumerate(w) if x==a] for a in set(w)}
 return [(a,b) for a in sorted(p) for b in sorted(p) if a<b and ((p[a][0]<p[b][0]<p[a][1]<p[b][1]) or (p[b][0]<p[a][0]<p[b][1]<p[a][1]))]
def rank(rows):
 b={}
 for raw in rows:
  r={j:F(a) for j,a in enumerate(raw) if a}
  while r:
   p=min(r)
   if p not in b:
    d=r[p];b[p]={j:a/d for j,a in r.items()};break
   d=r[p]
   for j,a in b[p].items():
    r[j]=r.get(j,0)-d*a
    if not r[j]:del r[j]
 return len(b)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--build-dir',type=Path,default=ROOT/'build');ap.add_argument('--matrix',action='store_true');args=ap.parse_args();build=args.build_dir.resolve()
 words=load('results/mutants11/words.json');assert len({canonical(w) for w in words.values()})==4;assert all(graph(w)==graph(words['D']) for w in words.values())
 ordinary=[poly('results/mutants11/'+s+'_ordinary.json') for s in words];primitive=[poly('results/mutants11/'+s+'_primitive.json') for s in words];delta=poly('results/mutants11/difference.json')
 assert len(delta)==27 and all(len(p)==521 for p in ordinary) and all(len(p)==138 for p in primitive)
 assert sub(ordinary[0],ordinary[1])==delta==sub(primitive[0],primitive[1]);assert ordinary[0]==ordinary[2] and ordinary[1]==ordinary[3];assert primitive[0]==primitive[2] and primitive[1]==primitive[3]
 assert all(sub(a,b)==sub(ordinary[0],primitive[0]) for a,b in zip(ordinary,primitive))
 def at(p):return sum(c for k,c in p.items() if not any(k[j] for j in range(1,len(k)) if j!=2))
 assert [at(p) for p in ordinary]==[11342,11322,11342,11322];assert [at(p) for p in primitive]==[-4888,-4908,-4888,-4908]
 report={'mutants11':True,'ordinary_terms':521,'primitive_terms':138,'difference_terms':27,'primitive_certificates':{}}
 all_words=[];expected=[]
 for n in [11,12,13]:
  d=load(f'certificates/order{n}/certificate.json');es=d['entries'];alpha=[F(e['coefficient']) for e in es]
  hs=[list(map(int,e['sl2_primitive_scaled'] if n==11 else e['primitive_scaled'])) for e in es];m=(n-1)//2;fs=[F(h[m] if len(h)>m else 0,8**(n-m)) for h in hs];target=sum(a*b for a,b in zip(alpha,fs));assert target==F(d['target_sum'])!=0
  if n in (11,13):
   names=['K_Q_monomial','K_P_monomial','K_S_binomial'] if n==11 else ['Q','P','S'];sizes=[n+1,n+1,n+2]
   features=[sum(([int(a) for a in e[name]]+[0]*(size-len(e[name])) for name,size in zip(names,sizes)),[]) for e in es]
   r=17 if n==11 else 22;assert rank(features)==r and rank([row+[f] for row,f in zip(features,fs)])==r+1
   assert all(sum(a*row[j] for a,row in zip(alpha,features))==0 for j in range(sum(sizes)))
  report['primitive_certificates'][str(n)]={'diagrams':len(es),'exact_value':str(target)}
  all_words += [e['word'] for e in es];expected += hs
  if n==13:
   leaf14=[{'word':[13,e['word'][0],13]+e['word'][1:],'h':[-4*a for a in h]} for e,h in zip(es,hs)]
   for x in leaf14:all_words.append(x['word']);expected.append(x['h'])
 proc=subprocess.run([str(build/'sl2_word_batch')],input='\n'.join(' '.join(map(str,w)) for w in all_words)+'\n',text=True,capture_output=True,check=True)
 fresh=[json.loads(s) for s in proc.stdout.splitlines()];assert len(fresh)==len(expected)
 assert all(list(map(int,d['primitive_scaled']))==h for d,h in zip(fresh,expected));report['fresh_sl2_polynomials']=len(fresh)
 d=load('certificates/order13/certificate.json');es=d['entries'];sp=[json.loads(s) for s in (ROOT/'certificates/order13/certificate_spin.ndjson').read_text().splitlines()]
 for r,e in zip(sp,es):
  assert all(sum(a*math.comb(N,k) for k,a in enumerate(e['S']) if k<=N)==int(v) for N,v in enumerate(r['primitive_S_values_N0_to16']))
 if args.matrix:
  for a in range(1,9):
   r=subprocess.run([str(build/'sl2_direct_matrix'),str(ROOT/'certificates/order13/certificate_input.json'),str(a)],text=True,capture_output=True,check=True);rows=[json.loads(s) for s in r.stdout.splitlines()];assert len(rows)==23 and all(x['passed'] for x in rows)
  report['independent_matrix_checks']=184
 k=poly('results/K11/K11.json');assert len(k)==587;report['K11_terms']=587
 for name,m in [('standard',lambda e:(sum((j-1)*e[j] for j in range(1,len(e)))+e[0],)),('scalar',lambda e:(e[0]+sum(e[1:]),sum(j*e[j] for j in range(1,len(e)))) )]:
  out={}
  for e,c in k.items():q=m(e);out[q]=out.get(q,0)+c
  out={q:c for q,c in out.items() if c};assert out==({(1,):1} if name=='standard' else {(11,22):1})
 cert=load('results/mutants14/certificate.json');assert len(set(cert['exact_integer_weights']))==4;assert len({canonical(w) for w in cert['words']})==4;assert cert['CRT_product_greater_than_twice_bound'];assert all(graph(w)==graph(cert['words'][0]) for w in cert['words'])
 report['mutants14']=True;report['passed']=True
 out=ROOT/'verification.json';out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');print(json.dumps(report,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
