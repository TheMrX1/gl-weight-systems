import ast
from collections import defaultdict
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time
import networkx as nx
from check_chromatic import DATA, four_terms, graph

HERE=Path(__file__).resolve().parent


def polynomial(text):
                                                                         
    def add(a,b,sign=1):
        p=a.copy()
        for k,v in b.items(): p[k]=p.get(k,Fraction(0))+sign*v
        return {k:v for k,v in p.items() if v}
    def visit(node):
        if isinstance(node,ast.Constant) and type(node.value)==int:
            return {0:Fraction(node.value)}
        if isinstance(node,ast.Name) and node.id=='c':return {1:Fraction(1)}
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):
            return {k:-v for k,v in visit(node.operand).items()}
        if isinstance(node,ast.BinOp):
            a,b=visit(node.left),visit(node.right)
            if isinstance(node.op,ast.Add):return add(a,b)
            if isinstance(node.op,ast.Sub):return add(a,b,-1)
            if isinstance(node.op,ast.Mult):
                p=defaultdict(Fraction)
                for i,u in a.items():
                    for j,v in b.items():p[i+j]+=u*v
                return {k:v for k,v in p.items() if v}
            if isinstance(node.op,ast.Div) and set(b)=={0}:
                return {k:v/b[0] for k,v in a.items()}
            if isinstance(node.op,ast.Pow) and isinstance(node.right,ast.Constant):
                exponent=node.right.value
                if not isinstance(exponent,int) or not 0<=exponent<=20:raise ValueError('Bad exponent')
                if a=={1:Fraction(1)}:return {exponent:Fraction(1)}
        raise ValueError('Unsupported polynomial syntax')
    return visit(ast.parse(text.replace('^','**'),mode='eval').body)


def main():
    start=time.monotonic()
    with (DATA/'cert.9.int.csv').open() as f: circles=list(csv.DictReader(f))
    with (DATA/'cert.9.csv').open() as f: relations=list(csv.DictReader(f))
    graphs=[];buckets=defaultdict(list);lookup={};words=[]
    def find(text):
        if text in lookup:return lookup[text]
        g=graph(text)
                                                                                 
        key=nx.weisfeiler_lehman_graph_hash(g,iterations=3)
        for i in buckets[key]:
            if nx.is_isomorphic(g,graphs[i]):
                lookup[text]=i
                return i
        i=len(graphs);graphs.append(g);buckets[key].append(i);lookup[text]=i
        return i
    for row in circles:
        i=find(row['graph6'])
        assert i==len(words),'Duplicate author intersection graph'
        words.append(row)
    coefficients=defaultdict(Fraction)
    traces=[]
    for j,row in enumerate(relations):
        ids=[]
        for sign,text in zip([1,-1,1,-1],four_terms(graph(row['graph6']),int(row['x']),int(row['y']))):
            i=find(text);coefficients[i]+=sign*Fraction(row['coef']);ids.append(i)
        traces.append(ids)
        if (j+1)%500==0:print('4T',j+1,'isomorphism classes',len(graphs),flush=True)
    nonzero={i:c for i,c in coefficients.items() if c}
    assert len(nonzero)==len(words)==5006
    assert all(i<len(words) for i in nonzero),'A non-circle graph survived'
                                                                         
    total=defaultdict(Fraction)
    for i,c in nonzero.items():
        for k,v in polynomial(words[i]['sl2']).items():total[k]+=c*v
    total={k:v for k,v in total.items() if v}
    expected={4:Fraction(1),3:Fraction(-41,32),2:Fraction(63,256),1:Fraction(9,256)}
    assert total==expected,(total,expected)
    entries=[{'id':i,'coefficient':str(nonzero[i]),'graph6':row['graph6'],
              'word':list(map(int,row['chord diagram'].split()))} for i,row in enumerate(words)]
    report={'repository_commit':json.loads((DATA/'SOURCE.json').read_text())['commit'],
            'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [DATA/'cert.9.csv',DATA/'cert.9.int.csv']},
            'relations':len(relations),'terms_before_isomorphism':4*len(relations),
            'isomorphism_classes_encountered':len(graphs),'surviving_circle_graphs':len(entries),
            'surviving_non_circle_graphs':0,
            'author_sl2_polynomial_sum':{str(k):str(v) for k,v in sorted(total.items())},
            'author_sl2_polynomials_recomputed':False,
            'elapsed_seconds':time.monotonic()-start,'entries':entries,'relation_term_ids':traces}
    (HERE.parent/'certificates/order9/reconstructed.json').write_text(json.dumps(report,indent=2)+'\n')
    print({k:v for k,v in report.items() if k not in ['entries','relation_term_ids']})


if __name__=='__main__':main()
