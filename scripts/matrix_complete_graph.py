import json
import subprocess
from pathlib import Path
import numpy as np

LIMIT=(1<<63)-1

def dot(a,b):
    bound=int(np.max(np.abs(a)))*int(np.max(np.abs(b)))*a.shape[1]
    if a.dtype==object or b.dtype==object or bound>LIMIT:
        return a.astype(object) @ b.astype(object)
    return a @ b

def character(N,t,max_casimir):
    d=N*N
    units=[]
    for i in range(N):
        row=[]
        for j in range(N):
            x=np.zeros((d,d),dtype=np.int64)
            for k in range(N):
                for l in range(N):
                    if j==k:x[i*N+l,k*N+l]+=1
                    if l==i:x[k*N+j,k*N+l]-=1
            if i==j:x+=t*np.eye(d,dtype=np.int64)
            row.append(x)
        units.append(row)
    B=np.block(units)
    b=np.eye(N*d,dtype=np.int64)
    casimirs=[N]
    for k in range(1,max_casimir+1):
        b=dot(b,B)
        trace=sum(int(x) for x in b.diagonal())
        value,remainder=divmod(trace-N*t**k,d-1)
        assert remainder==0
        casimirs.append(value)
    omega=np.zeros((d*d,d*d),dtype=np.int64)
    for i in range(N):
        for j in range(N):omega+=np.kron(units[i][j],units[j][i])
    return casimirs,omega

def weight(n,N,t,max_casimir):
    casimirs,omega=character(N,t,max_casimir)
    d=N*N
    power=np.eye(d*d,dtype=np.int64)
    for _ in range(n):power=dot(power,omega)
    trace=sum(int(power[b*d+a,a*d+b]) for a in range(d) for b in range(d))
    value,remainder=divmod(trace-(N*t*t)**n,d-1)
    assert remainder==0
    return casimirs,value

def evaluate(terms,casimirs):
    total=0
    for term in terms:
        value=int(term['coefficient'])
        for x,e in zip(casimirs,term['powers']):value*=x**e
        assert all(e==0 for e in term['powers'][len(casimirs):])
        total+=value
    return total

def selftest():
    exe=Path(__file__).resolve().parents[1]/'build/symbolic_witness'
    for n in range(1,6):
        word=' '.join(map(str,list(range(n))*2))+'\n'
        terms=json.loads(subprocess.run([str(exe),'--word'],input=word,capture_output=True,text=True,check=True).stdout)
        maxk=max(k for term in terms for k,e in enumerate(term['powers']) if e)
        for N,t in [(2,0),(3,1)]:
            cs,v=weight(n,N,t,maxk)
            assert evaluate(terms,cs)==v,(n,N,t)
        print('Matrix formula PASS K'+str(n),flush=True)

if __name__=='__main__':selftest()
