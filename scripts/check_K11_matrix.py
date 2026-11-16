from pathlib import Path
import json
from matrix_complete_graph import weight,evaluate
ROOT=Path(__file__).resolve().parents[1]
terms=json.loads((ROOT/'results/K11/K11.json').read_text())
checks=[]
for N in [2,3,4]:
 for t in [0,1,-1,2]:
  cs,value=weight(11,N,t,22)
  assert evaluate(terms,cs)==value
  checks.append({'N':N,'twist':t,'value':str(value)})
(ROOT/'results/K11/fresh_matrix_checks.json').write_text(json.dumps({'passed':True,'checks':checks},indent=2)+'\n')
print('PASS: 12 exact matrix contractions agree with K11')
