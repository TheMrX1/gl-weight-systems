from collections import defaultdict
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]/'results/hooks'


def add(a,b,sign=1):
    out=defaultdict(int,a)
    for p,c in b.items():out[p]+=sign*c
    return {p:c for p,c in out.items() if c}


def multiply(a,b):
    out=defaultdict(int)
    for p,c in a.items():
        for q,d in b.items():out[tuple(x+y for x,y in zip(p,q))]+=c*d
    return {p:c for p,c in out.items() if c}


def power(a,n):
    out={(0,0,0,0,0):1}
    for _ in range(n):out=multiply(out,a)
    return out


def scale(a,n):return {p:c*n for p,c in a.items() if c*n}


def main():
    names=['N','A','B','X','Y'];variables=[]
    for i in range(5):
        p=[0]*5;p[i]=1;variables.append({tuple(p):1})
    N,A,B,X,Y=variables
    casimirs=[N]+[add(multiply(A,power(X,m-1)),multiply(B,power(Y,m-1))) for m in range(1,7)]
    terms=json.loads((HERE/'universal_delta.json').read_text())
    original={}
    for term in terms:
        value={(0,0,0,0,0):int(term['coefficient'])}
        for m,e in enumerate(term['powers']):
            if e:value=multiply(value,power(casimirs[m],e))
        original=add(original,value)
    factor1=power(add(X,Y,-1),2)
    factor2=add(add(multiply(A,Y),multiply(B,X)),multiply(N,multiply(X,Y)),-1)
    factor3=add(add(scale(add(A,B),3),power(N,2)),multiply(N,add(X,Y)),-1)
    factored=scale(multiply(multiply(A,B),multiply(factor1,multiply(factor2,power(factor3,2)))),20)
    assert original==factored
    certificate={'variables':names,'identity':'Delta=20*A*B*(X-Y)^2*(A*Y+B*X-N*X*Y)*(3*A+3*B+N^2-N*X-N*Y)^2',
                 'expanded_term_count':len(original),'all_coefficients_identical':True,
                 'uses_sympy':False,
                 'expanded_terms':[{'powers':list(p),'coefficient':str(c)} for p,c in sorted(original.items())],
                 'universal_delta_sha256':hashlib.sha256((HERE/'universal_delta.json').read_bytes()).hexdigest()}
    (HERE/'two_atom_identity_certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    print('PASS: independent sparse multiplication,',len(original),'identical coefficients')


if __name__=='__main__':main()
