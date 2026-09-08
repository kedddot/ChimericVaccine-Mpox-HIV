"""
AllerTOP v2 reimplementation probe.

Method as published (Dimitrov, Bangov, Flower & Doytchinova 2014, J Mol Model
20:2278):
  1. each residue -> 5 E-descriptors (Venkatarajan & Braun 2001)
  2. auto- and cross-covariance (ACC) over lags 1..L  -> 25*L features
  3. k-nearest-neighbour classification against 2,427 allergens +
     2,427 non-allergens
  Published 5-fold cross-validation accuracy: ~88.7%

We do NOT know their exact lag, k, or preprocessing. So this script sweeps them
and asks which combination reproduces the published accuracy. Hitting ~88-89%
is evidence the descriptor table and ACC formulation are right; missing it
badly means they are not.
"""
import numpy as np, itertools, sys

# Venkatarajan & Braun (2001) E1-E5.
E = {
'A':(0.008,0.134,-0.475,-0.039,0.181),  'R':(0.171,-0.361,0.107,-0.258,-0.364),
'N':(0.255,0.038,0.117,0.118,-0.055),   'D':(0.303,-0.057,-0.014,0.225,0.156),
'C':(-0.132,0.174,0.070,0.565,-0.374),  'Q':(0.149,-0.184,-0.030,0.035,-0.112),
'E':(0.221,-0.280,-0.315,0.157,0.303),  'G':(0.218,0.562,-0.024,0.018,0.106),
'H':(0.023,-0.177,0.041,0.280,-0.021),  'I':(-0.353,0.071,-0.088,-0.195,-0.107),
'L':(-0.267,0.018,-0.265,-0.274,0.206), 'K':(0.243,-0.339,-0.044,-0.325,-0.027),
'M':(-0.239,-0.141,-0.155,0.321,0.077), 'F':(-0.329,-0.023,0.072,-0.002,0.208),
'P':(0.173,0.286,0.407,-0.215,0.384),   'S':(0.199,0.238,-0.015,-0.068,-0.196),
'T':(0.068,0.147,-0.015,-0.132,-0.274), 'W':(-0.296,-0.186,0.389,0.083,0.297),
'Y':(-0.141,-0.057,0.425,-0.096,-0.091),'V':(-0.274,0.136,-0.187,-0.196,-0.299),
}

def readfa(p):
    recs, h, cur = [], None, []
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if h is not None: recs.append((h, "".join(cur)))
            h, cur = line[1:], []
        else:
            cur.append(line.strip())
    if h is not None: recs.append((h, "".join(cur)))
    return recs

def acc(seq, L):
    """ACC over lags 1..L -> 25*L vector."""
    m = np.array([E[c] for c in seq if c in E], dtype=float)   # N x 5
    n = len(m)
    if n == 0:
        return np.zeros(25 * L)
    out = []
    for lag in range(1, L + 1):
        if n - lag < 1:
            out.append(np.zeros(25))
            continue
        a = m[:n - lag]          # (n-lag) x 5
        b = m[lag:]              # (n-lag) x 5
        # Z[j,k] = sum_i E_j(i) * E_k(i+lag) / (n-lag)
        z = (a.T @ b) / (n - lag)
        out.append(z.ravel())
    return np.concatenate(out)

def featurize(seqs, L):
    return np.vstack([acc(s, L) for s in seqs])

def cv_accuracy(X, y, k, folds=5, seed=0):
    rng = np.random.RandomState(seed)
    idx = rng.permutation(len(y))
    accs = []
    for f in range(folds):
        te = idx[f::folds]
        tr = np.setdefault = np.array([i for i in idx if i not in set(te.tolist())])
        Xtr, ytr, Xte, yte = X[tr], y[tr], X[te], y[te]
        # euclidean distances
        d = ((Xte**2).sum(1)[:, None] + (Xtr**2).sum(1)[None, :]
             - 2 * Xte @ Xtr.T)
        nn = np.argpartition(d, k, axis=1)[:, :k]
        pred = (ytr[nn].mean(1) > 0.5).astype(int)
        accs.append((pred == yte).mean())
    return float(np.mean(accs))

if __name__ == "__main__":
    al = readfa("/Users/nek/Downloads/allergens_dataset.fasta")
    na = readfa("/Users/nek/Downloads/non-allergens_dataset.fasta")
    seqs = [s for _, s in al] + [s for _, s in na]
    y = np.array([1]*len(al) + [0]*len(na))
    print(f"training set: {len(al)} allergens + {len(na)} non-allergens")

    print(f"\n{'lag':>4} {'feats':>6} " + " ".join(f"k={k:<2d}" for k in (1,3,5)))
    for L in (1, 3, 5, 7, 9):
        X = featurize(seqs, L)
        # AllerTOP scales descriptors; test both raw and standardised
        row = f"{L:>4} {X.shape[1]:>6} "
        for k in (1, 3, 5):
            row += f"{100*cv_accuracy(X, y, k):5.1f} "
        print(row + "   (raw)")
        Xs = (X - X.mean(0)) / (X.std(0) + 1e-12)
        row = f"{'':>4} {'':>6} "
        for k in (1, 3, 5):
            row += f"{100*cv_accuracy(Xs, y, k):5.1f} "
        print(row + "   (standardised)")
