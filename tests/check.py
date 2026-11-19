#!/usr/bin/env python3
"""Check a bandLD query result against the exact reference computed from the generator's dosage numerators.
usage: check.py DATA_PREFIX MATRIX_FILE WINDOW_BP ROUND_BITS
DATA_PREFIX.k.npy / .miss.npy / .z come from tests/make_bgen.py; MATRIX_FILE.snps.tsv says which SNPs are in the
matrix and in which allele orientation (the generator's alleles are A/G, dosage of G). The matrix must have been
written with --precision 17. Exact mode: every stored value must be the same double as the reference."""
import gzip, math, sys
import numpy as np
pre, mat, W, B = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
K = np.load(pre + ".k.npy").astype(object)          # Python integers: exact
N, M = K.shape
z = [l.split() for l in open(pre + ".z").read().strip().split("\n")[1:]]
ids = [r[0] for r in z]; pos = {r[0]: int(r[2]) for r in z}; col = {r: i for i, r in enumerate(ids)}
snps = [l.rstrip("\n").split("\t") for l in open(mat + ".snps.tsv") if not l.startswith("#")]
op = gzip.open if mat.endswith(".gz") else open
rows = [l.split() for l in op(mat, "rt").read().strip().split("\n")]
n = len(snps)
assert len(rows) == n and all(len(r) == n for r in rows), "matrix shape"
sx = [int(K[:, j].sum()) for j in range(M)]; sxx = [int((K[:, j] * K[:, j]).sum()) for j in range(M)]
V = [N * sxx[j] - sx[j] * sx[j] for j in range(M)]
bad = 0; nstored = 0; nnan = 0
for x in range(n):
    jx = col[snps[x][2]]; fx = -1 if snps[x][4] != "G" else 1
    assert rows[x][x] == "1", "diagonal"
    for y in range(x):
        jy = col[snps[y][2]]; fy = -1 if snps[y][4] != "G" else 1
        got = rows[x][y]
        assert rows[y][x] == got, "symmetry"
        if abs(pos[snps[x][2]] - pos[snps[y][2]]) > W:
            nnan += 1
            if got != "nan": bad += 1
            continue
        nstored += 1
        C = int((K[:, jx] * K[:, jy]).sum())
        r = float(N * C - sx[jx] * sx[jy]) / math.sqrt(float(V[jx]) * float(V[jy]))
        if B: r = round(r * 2**B) / 2**B if abs(r * 2**B - round(r * 2**B)) != 0.5 else float(got) * fx * fy   # ties: either way
        r *= fx * fy
        if r == 0: r = 0.0
        if float(got) != r:
            bad += 1
            if bad <= 5: print("MISMATCH", snps[x][2], snps[y][2], got, repr(r))
print(f"{mat}: {n} SNPs, {nstored} pairs within the window checked, {nnan} beyond it (nan), mismatches {bad}")
sys.exit(1 if bad else 0)
