import json
import math
import subprocess
import time
import os
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/K11'
EXE=ROOT/'build/symbolic_traceless_compressed'
OLD=ROOT/'build/symbolic_witness'
ZERO=(0,)*23

def add(out,p,scale=1):
    for k,c in p.items():
        c=out.get(k,0)+scale*c
        if c:out[k]=c
        else:out.pop(k,None)

def mul(a,b):
    out={}
    for k,c in a.items():
        for l,d in b.items():
            m=tuple(x+y for x,y in zip(k,l))
            out[m]=out.get(m,0)+c*d
    return {k:c for k,c in out.items() if c}

def shifted(j):
    out={}
    for l in range(j+1):
        key=[0]*23
        key[0]=-(j-l)
        key[1]=j-l
        key[l]+=1
        add(out,{tuple(key):math.comb(j,l)*(-1)**(j-l)})
    return out

def parse(terms):
    out={tuple(t['powers']):int(t['coefficient']) for t in terms}
    assert len(out)==len(terms)
    return {k:c for k,c in out.items() if c}

def serialize(p):
    return [{'coefficient':str(c),'powers':list(k)} for k,c in sorted(p.items(),key=lambda pair:pair[0][::-1])]

def reconstruct(fs,n):
    out={}
    powers={(j,0):{ZERO:1} for j in range(2,23)}
    def power(j,e):
        if (j,e) not in powers:powers[j,e]=mul(power(j,e-1),shifted(j))
        return powers[j,e]
    for k,f in enumerate(fs[:n+1]):
        for m,c in f.items():
            assert m[1]==0
            base=[0]*23;base[0]=m[0]-(n-k);base[1]=2*(n-k)
            p={tuple(base):c*math.comb(n,k)}
            for j,e in enumerate(m[2:],2):
                if e:p=mul(p,power(j,e))
            add(out,p)
    assert all(min(k)>=0 for k in out),'negative Laurent powers did not cancel'
    return out

def save(path,obj):
    if path.exists():
        assert json.loads(path.read_text())==obj,('refusing to overwrite different result',path)
        return
    with path.open('x') as f:json.dump(obj,f);f.write('\n')

def main():
    start=time.monotonic()
    OUT.mkdir(exist_ok=True)
    fs=[{ZERO:1}]
    timings=[]
    for n in range(1,12):
        path=OUT/f'F{n}.json'
        text=' '.join(map(str,list(range(n))*2))+'\n'
        t=time.monotonic()
        if path.exists():
            f=parse(json.loads(path.read_text()))
        else:
            print('Starting C1=0 K'+str(n),flush=True)
            with (OUT/f'F{n}_progress.log').open('a') as log:
                env=dict(os.environ)
                if n==11:
                    env['CHORD_DAG_CHECKPOINT']=str(OUT/'F11.dag')
                    if not (OUT/'F11.dag').exists():
                        build_env=dict(env,CHORD_DAG_BUILD_ONLY='1')
                        subprocess.run([str(EXE),'--word'],input=text,text=True,
                                       stdout=subprocess.DEVNULL,stderr=log,check=True,timeout=3600,env=build_env)
                run=subprocess.run([str(EXE),'--word'],input=text,text=True,
                                   stdout=subprocess.PIPE,stderr=log,check=True,timeout=3600,env=env)
            f=parse(json.loads(run.stdout))
            save(path,serialize(f))
        fs.append(f)
        assert all(k[1]==0 for k in f)
        p=reconstruct(fs,n)
        if n<=9:
            baseline=subprocess.run([str(OLD),'--word'],input=text,text=True,
                                    stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,check=True,timeout=180)
            assert p==parse(json.loads(baseline.stdout)),n
            print('Exact full-gl reconstruction PASS K'+str(n),flush=True)
        save(OUT/f'K{n}.json',serialize(p))
        duration=time.monotonic()-t
        timings.append({'n':n,'seconds':duration,'traceless_terms':len(f),'full_terms':len(p)})
        print('Completed K'+str(n),timings[-1],flush=True)
    if not (OUT/'timings.json').exists():
        save(OUT/'timings.json',{'runs':timings,'total_seconds':time.monotonic()-start})

if __name__=='__main__':main()
