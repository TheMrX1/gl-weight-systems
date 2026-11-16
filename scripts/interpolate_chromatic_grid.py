import json
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]/'results/graphs9'
P=1000000007

def interpolate(points):
    result=[0]*len(points)
    for i,(x,y) in enumerate(points):
        basis=[1]
        denom=1
        for j,(xx,_) in enumerate(points):
            if i==j: continue
            new=[0]*(len(basis)+1)
            for k,a in enumerate(basis):
                new[k]=(new[k]-xx*a)%P
                new[k+1]=(new[k+1]+a)%P
            basis=new
            denom=denom*(x-xx)%P
        scale=y*pow(denom,-1,P)%P
        result=[(a+scale*b)%P for a,b in zip(result,basis)]
    return result

def evaluate(poly,x):
    r=0
    for a in reversed(poly):r=(r*x+a)%P
    return r

rows=[json.loads(s) for s in (HERE/'grid_values.ndjson').read_text().splitlines()]
assert rows[-1].get('summary') and len(rows)==5007
entries=json.loads((HERE.parents[1]/'certificates/order9/certificate.json').read_text())['entries']
assert len(entries)==5006
totals=[0]*110
for row,entry in zip(rows[:-1],entries):
    assert row['id']==entry['id']
    alpha=Fraction(entry['coefficient'])
    alpha_mod=alpha.numerator*pow(alpha.denominator,-1,P)%P
    for i,value in enumerate(row['values']):totals[i]=(totals[i]+alpha_mod*value)%P
lanes=rows[-1]['lanes']
assert totals==[lane['value'] for lane in lanes]
lookup={(lane['N'],lane['x']):value for lane,value in zip(lanes,totals)}
by_n={N:interpolate([(x,lookup[N,x]) for x in range(19)]) for N in range(1,6)}
for N,poly in by_n.items():assert evaluate(poly,19)==lookup[N,19]
by_g=[[0]*19 for _ in range(5)]
for x_degree in range(19):
    n_poly=interpolate([(N*N,poly[x_degree]*pow(N,-1,P)%P) for N,poly in by_n.items()])
    for j,coef in enumerate(n_poly):by_g[4-j][x_degree]=coef
for x in range(1,11):
    actual=sum(evaluate(poly,x)*pow(6,9-2*g,P) for g,poly in enumerate(by_g))%P
    assert actual==lookup[6,x]
assert not any(by_g[0]), 'Chromatic graph 4T control failed'
assert all(not any(poly[10:]) for poly in by_g), 'Unexpected high x-degree'
expected={1:810,2:15930,3:-225360,4:594720}
for g,value in expected.items():assert evaluate(by_g[g],4)==value%P

terms=[]
for g,poly in enumerate(by_g):
    for k,a in enumerate(poly):
        if a:terms.append({'genus':g,'x_degree':k,'coefficient_mod_p':a,
                           'symmetric_representative':a if a<=P//2 else a-P})
report={'prime':P,'normalization':'X=N^(-9)*W(C_j=x*N^(j-1))',
        'entries':5006,'nonzero_coefficients_mod_p':terms,
        'coefficient_status':'Nonzero residues are rigorous nonvanishing witnesses; zero residues alone do not prove exact vanishing.',
        'interpolation_bounds':'degree_x<=18 by PBW filtration; W odd degree_N<=9 by prechromatic theorem.',
        'checks':{'independent_sum':True,'five_x19_holdouts':True,
                  'ten_N6_holdouts':True,'chromatic_part_zero_mod_p':True,
                  'old_exact_x4_coefficients':True}}
(HERE/'bivariate_certificate_mod_p.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False,indent=2))
