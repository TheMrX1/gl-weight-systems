from pathlib import Path
import argparse,json,subprocess,math
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--recompute',action='store_true');args=ap.parse_args()
cert=ROOT/'certificates/order13';es=json.loads((cert/'certificate.json').read_text())['entries']
if args.recompute:
 outputs={}
 for exe in ['graph_corank','graph_diagonal','graph_pairs']:
  p=subprocess.run([str(ROOT/'build'/exe),str(cert/'certificate_input.json')],capture_output=True,text=True,check=True)
  outputs[exe]=[json.loads(s) for s in p.stdout.splitlines()]
 for e,q,p,s in zip(es,outputs['graph_corank'],outputs['graph_diagonal'],outputs['graph_pairs']):
  assert q['coefficients']==e['Q'] and p['coefficients']==e['P']
  assert all(sum(a*math.comb(N,k) for k,a in enumerate(e['S']) if k<=N)==int(v) for N,v in enumerate(s['primitive_S_values_N0_to16']))
else:
 for filename,key,family in [('certificate_corank.ndjson','corank_primitive','Q'),('certificate_orthogonal.ndjson','orthogonal_primitive','P')]:
  rows={r['id']:r[key] for line in (cert/filename).read_text().splitlines() if key in (r:=json.loads(line))}
  assert len(rows)==23 and all(rows[i]==e[family] for i,e in enumerate(es))
 rows=[json.loads(line) for line in (cert/'certificate_spin.ndjson').read_text().splitlines()]
 assert len(rows)==23
 for i,(r,e) in enumerate(zip(rows,es)):
  assert r['id']==i
  assert all(sum(a*math.comb(N,k) for k,a in enumerate(e['S']) if k<=N)==int(v) for N,v in enumerate(r['primitive_S_values_N0_to16']))
print('PASS: 23 full Q/P/S polynomial triples', 'recomputed' if args.recompute else 'checked from saved data')
