import csv
import itertools
import json
import random
import subprocess
import time
from fractions import Fraction
from pathlib import Path

import networkx as nx

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / 'certificates/order9/source'


def graph(text):
    return nx.from_graph6_bytes(text.encode())


def encode(g):
    return nx.to_graph6_bytes(g, header=False).decode().strip()


def four_terms(g, x, y):
    switched = g.copy()
    for v in g.neighbors(x):
        if v != y:
            if switched.has_edge(y, v):
                switched.remove_edge(y, v)
            else:
                switched.add_edge(y, v)
    both = switched.copy()
    edge = g.copy()
    for h in [both, edge]:
        if h.has_edge(x, y): h.remove_edge(x, y)
        else: h.add_edge(x, y)
    return [encode(h) for h in [g, switched, both, edge]]


def value(p, q):
    return sum(c*q**k for k, c in enumerate(p))


def direct_count(g, q):
    adjacency = [list(g.neighbors(i)) for i in range(len(g))]
    colors = [-1]*len(g)
    def visit(i):
        if i == len(g): return 1
        total = 0
        for color in range(q):
            if all(colors[j] != color for j in adjacency[i] if j < i):
                colors[i] = color
                total += visit(i+1)
        colors[i] = -1
        return total
    return visit(0)


def partitions(n):
    def rec(i, blocks):
        if i == n:
            yield [b[:] for b in blocks]
            return
        for b in blocks:
            b.append(i)
            yield from rec(i+1, blocks)
            b.pop()
        blocks.append([i])
        yield from rec(i+1, blocks)
        blocks.pop()
    yield from rec(0, [])


def main():
    started = time.monotonic()
    with (DATA/'cert.9.csv').open() as f: rows = list(csv.DictReader(f))
    with (DATA/'cert.9.int.csv').open() as f: circles = list(csv.DictReader(f))
    assert len(rows) == 3300 and len(circles) == 5006
    relations = [four_terms(graph(r['graph6']), int(r['x']), int(r['y'])) for r in rows]
    requested = {s for rel in relations for s in rel}
    requested.update(r['graph6'] for r in circles)
    rng = random.Random(20261101)
    random_relations = []
    for _ in range(300):
        n = rng.randrange(2, 10)
        g = nx.Graph(); g.add_nodes_from(range(n))
        g.add_edges_from(e for e in itertools.combinations(range(n), 2) if rng.randrange(2))
        x, y = rng.sample(range(n), 2)
        random_relations.append(four_terms(g, x, y))
    requested.update(s for rel in random_relations for s in rel)
    small = []
    for n in range(6):
        edges = list(itertools.combinations(range(n), 2))
        for bits in range(1 << len(edges)):
            g = nx.Graph(); g.add_nodes_from(range(n))
            g.add_edges_from(e for k, e in enumerate(edges) if bits >> k & 1)
            small.append(encode(g))
    requested.update(small)
    source = '\n'.join(sorted(requested))+'\n'
    run = subprocess.run([str(HERE.parent/'build/chromatic')],input=source,text=True,capture_output=True,check=True)
    polys = {}
    for line in run.stdout.splitlines():
        s, *terms = line.split(); polys[s] = list(map(int, terms))
    assert len(polys) == len(requested)
    combined = [Fraction(0)]*10
    for r, rel in zip(rows, relations):
        difference = [sum(sign*polys[s][k] for sign, s in zip([1,-1,1,-1], rel)) for k in range(10)]
        assert not any(difference)
        for k in range(10): combined[k] += Fraction(r['coef'])*difference[k]
    assert not any(combined)
    for rel in random_relations:
        for k in range(len(polys[rel[0]])):
            assert sum(sign*polys[s][k] for sign,s in zip([1,-1,1,-1],rel)) == 0
    checks = 0
    for s in small:
        g = graph(s)
        for q in range(4):
            assert value(polys[s], q) == direct_count(g, q)
            checks += 1
                                                                    
                                                                      
    primitive_checks = 0
    from math import factorial
    for s in small[::31]:
        g = graph(s)
        if not len(g): continue
        for q in [2,3]:
            total = 0
            for blocks in partitions(len(g)):
                term = (-1)**(len(blocks)-1)*factorial(len(blocks)-1)
                for block in blocks:
                    h = nx.convert_node_labels_to_integers(g.subgraph(block))
                    term *= direct_count(h, q)
                total += term
            assert total == q*polys[s][1]
            primitive_checks += 1
                                                                         
    for row in circles:
        word = list(map(int,row['chord diagram'].split()))
        assert len(word)==18 and set(word)==set(range(9))
        pos = {k:[i for i,x in enumerate(word) if x==k] for k in range(9)}
        assert all(len(p)==2 for p in pos.values())
        h = nx.Graph(); h.add_nodes_from(range(9))
        for a,b in itertools.combinations(range(9),2):
            i,j=pos[a]; k,l=pos[b]
            if (i<k<j) != (i<l<j): h.add_edge(a,b)
        assert nx.is_isomorphic(graph(row['graph6']),h)
    output = {
        'source_commit':json.loads((DATA/'SOURCE.json').read_text())['commit'],
        'four_term_relations':len(rows), 'supplied_intersection_graphs_checked':len(circles),
        'labeled_graph_polynomials_computed':len(polys),
        'chromatic_combination_coefficients_q0_to_q9':list(map(str,combined)),
        'primitive_chromatic_combination_coefficients_q0_to_q9':['0']*10,
        'primitive_combination_method':'chi(pi(G)) = q*[q]chi(G); applied linearly to the verified 4T sum',
        'exhaustive_labeled_graphs_through_order5':len(small),
        'independent_coloring_count_checks':checks,
        'random_four_term_checks':len(random_relations),
        'independent_primitive_partition_checks':primitive_checks,
        'limitations':'Does not recompute author sl2 polynomials, canonical 5006-term expansion, or rectangular gl weights on the certificate.',
        'elapsed_seconds':time.monotonic()-started,
    }
    (HERE/'chromatic_report.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__': main()
