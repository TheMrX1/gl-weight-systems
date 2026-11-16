import json
from pathlib import Path
from check_sparse_identity import add,multiply,power,scale

HERE=Path(__file__).resolve().parents[1]/'results/hooks'
ONE={(0,0,0,0,0):1}


def main():
    variables=[]
    for i in range(3):
        p=[0]*5;p[i]=1;variables.append({tuple(p):1})
    terms=json.loads((HERE/'universal_delta.json').read_text())
    results=[]
    for n in [1,2,3]:
                                                                                 
        rows=[]
        for i in range(n):
            row={}
            for j in range(i,n):row=add(row,variables[j])
            rows.append(row)
        series=[ONE]+[{} for _ in range(7)]
        for i,row in enumerate(rows):
            level=add(row,scale(ONE,n-i-1));numerator=add(level,ONE)
            old=series.copy()
            for k in range(1,8):
                series[k]=add(add(old[k],multiply(numerator,old[k-1]),-1),multiply(level,series[k-1]))
        c=[scale(ONE,n)]+[scale(series[m+1],-1) for m in range(1,7)]
        assert series[1]==scale(ONE,-n)
        result={}
        for t in terms:
            value=scale(ONE,int(t['coefficient']))
            for k,e in enumerate(t['powers']):
                if e:value=multiply(value,power(c[k],e))
            result=add(result,value)
        assert result=={}
        results.append({'rank':n,'highest_weights':'lambda_i=x_i+...+x_N',
                        'HC_difference_identically_zero':True,
                        'casimirs':[[{'powers':list(p[:n]),'coefficient':str(co)} for p,co in sorted(poly.items())] for poly in c]})
    report={'uses_sympy':False,'scope':'Fixed pair D,Dprime; all central characters, not just hooks',
            'conclusion':'HC is injective, so W_glN(D)=W_glN(Dprime) in Z(U(glN)) for N=1,2,3; same for primitive weights',
            'ranks':results}
    (HERE/'small_rank_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: Harish-Chandra difference is zero in ranks 1,2,3, by independent sparse arithmetic')


if __name__=='__main__':main()
