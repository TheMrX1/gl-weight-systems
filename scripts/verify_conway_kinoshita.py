from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1] / "results/kinoshita-terasaka-and-conway"
NAMES = ("conway", "kinoshita-terasaka")


def read_polynomial(path):
    terms = json.loads(path.read_text())
    result = {tuple(row["powers"]): int(row["coefficient"]) for row in terms}
    assert len(result) == len(terms)
    for powers, coefficient in result.items():
        assert len(powers) == 23 and coefficient != 0
        assert all(isinstance(k, int) and k >= 0 for k in powers)
        assert powers[0] + sum(i * k for i, k in enumerate(powers[1:], 1)) <= 22
    return result


def evaluate(polynomial, values, prime=None):
    total = 0
    for powers, coefficient in polynomial.items():
        term = coefficient
        for i, exponent in enumerate(powers):
            if exponent:
                term *= (pow(values[i], exponent, prime) if prime else values[i] ** exponent)
        total += term
    return total % prime if prime else total


def graph_corank(word):
    positions = [[i for i, x in enumerate(word) if x == a] for a in range(11)]
    adjacency = []
    for start, end in positions:
        row = sum(1 << b for b, (u, v) in enumerate(positions)
                  if (start < u < end) != (start < v < end))
        adjacency.append(row)
    basis = {}
    for row in adjacency:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = row
                break
            row ^= basis[pivot]
    return 11 - len(basis)


def check_standard(polynomial, corank):
    specialized = defaultdict(int)
    for powers, coefficient in polynomial.items():
        degree = powers[0] + sum((i - 1) * k for i, k in enumerate(powers[1:], 1))
        specialized[degree] += coefficient
    assert {k: a for k, a in specialized.items() if a} == {corank: 1}


def text_form(polynomial):
    lines = []
    for powers, coefficient in sorted(polynomial.items()):
        factors = []
        for k, exponent in enumerate(powers):
            if exponent:
                variable = "N" if k == 0 else f"C{k}"
                factors.append(variable + (f"^{exponent}" if exponent > 1 else ""))
        if abs(coefficient) != 1 or not factors:
            factors.insert(0, str(abs(coefficient)))
        lines.append(("- " if coefficient < 0 else "+ ") + "*".join(factors))
    if not lines:
        return "0\n"
    if lines[0].startswith("+ "):
        lines[0] = lines[0][2:]
    return "\n".join(lines) + "\n"


def main():
    exact = {}
    results = {}
    for name in NAMES:
        path = HERE / (name + ".json")
        exact[name] = polynomial = read_polynomial(path)
        word = list(map(int, (HERE / (name + ".word")).read_text().split()))
        assert Counter(word) == dict.fromkeys(range(11), 2)
        numeric_name = "conway" if name == "conway" else "mutant"
        prior = json.loads((HERE / (numeric_name + "_numeric.json")).read_text())
        modular = [evaluate(polynomial, xs, p) for xs, p in zip(prior["casimirs"], prior["primes"])]
        assert modular == prior["values"]
        matrix = {}
        for representation in ("hook", "square"):
            values = []
            for column in (0, 7):
                sample = json.loads((HERE / f"{name}-{representation}-{column}.json").read_text())
                assert sample["prime"] == 0 and sample["input_column"] == column
                value = int(sample["value"] if representation == "hook" else sample["scalar"])
                assert evaluate(polynomial, sample["casimirs"]) == value
                values.append(value)
            assert values[0] == values[1]
            matrix[representation] = values[0]
        corank = graph_corank(word)
        check_standard(polynomial, corank)
        results[name] = {
            "terms": len(polynomial),
            "max_casimir_index": max(i for powers in polynomial for i, k in enumerate(powers) if i and k),
            "sha256": sha256(path.read_bytes()).hexdigest(),
            "prime_moduli": prior["primes"],
            "modular_values": modular,
            "gl4_matrix_weights": matrix,
            "checked_columns": [0, 7],
            "corank_over_GF2": corank,
        }
    a, b = (exact[name] for name in NAMES)
    difference = {k: a.get(k, 0)-b.get(k, 0) for k in a.keys() | b.keys()
                  if a.get(k, 0) != b.get(k, 0)}
    assert a == b and len(a) == len(b) == 587
    assert not difference
    assert (HERE / "difference.txt").read_text().strip() == "0"
    assert not json.loads((HERE / "difference.json").read_text())
    for name, polynomial in exact.items():
        assert (HERE / (name + ".txt")).read_text() == text_form(polynomial)
    recorded = json.loads((HERE / "verification.json").read_text())
    assert recorded["diagrams"] == results
    assert recorded["coefficients_equal"] and recorded["difference_terms"] == len(difference)
    print("587 exact coefficients coincide; modular and matrix data agree.")


if __name__ == "__main__":
    main()
