This is a certificate proving that the sl(2)-weight system does not extend to all graphs with 9 vertices.
It consists of 3300 4T-relations and resolves into a combination of 5006 intersection graphs.

# `cert.9.int.csv`
In `cert.9.int.csv`, precomputed values of the sl(2)-weight system on relevant intersection graphs are provided.
The Killing form is the chosen form for the weight system.
The Chmutov-Varchenko chord deletion relation is used to compute the values.

CSV fields:
- `graph6`: the intersection graph in the graph6 format;
- `chord diagram`: an arbitrary chord diagram realising this intersection graph;
- `sl2`: the computed value of the sl(2)-weight system.

# `cert.9.csv`
In `cert.9.csv`, a combination of 4T-relations with non-zero sum and all non-zero graphs among intersection graphs is provided.
Graphs are stored together with pairs of vertices `(x, y)`.
Each relation is obtained the following way:
- term 1, sign +1: the original graph;
- term 2, sign -1: the original graph where for each neighbour `v` of `x` except `y`, edge `(y, v)` is switched;
- term 3, sign +1: the original graph where for each neighbour `v` of `x` except `y`, edge `(y, v)` is switched and edge `(x, y)` is also switched;
- term 4, sign -1: the original graph where edge `(x, y)` is switched.

CSV fields:
- `graph6`: the initial graph in the graph6 format;
- `x`, `y`: the vertices chosen for this relation;
- `coef`: the coefficient associated with this relation.

# `check.sage`
In `check.sage`, a SageMath program verifying the result is implemented.
First it verifies the result using all the data, then it verifies that `cert.9.int.csv` is correct.
The latter step is significantly longer and might take over an hour.

To run it, one needs to have installed:
1. SageMath.
2. Python library `tqdm` (used to track progress).
If not available system-wide, on a Unix-like system it can be installed locally with:

```bash
python -m venv env
source ./env/bin/activate
pip install tqdm
```

Then the script can be executed with the `sage check.sage` command.

- Daniil Fomichev
