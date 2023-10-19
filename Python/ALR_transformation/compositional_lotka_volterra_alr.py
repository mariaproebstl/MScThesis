import numpy as np
import sys


def choose_denom(P):
    """Pick a denominator for additive log-ratio transformation.
    """
    np.seterr(divide="ignore", invalid="ignore")
    log_change = None
    for p in P:
        s = p.sum(axis=1,keepdims=True)
        s[s==0] = 1
        deltas = np.log( (p/s)[1:] ) - np.log( (p/s)[:-1] )
        if log_change is None:
            log_change = deltas
        else:
            log_change = np.vstack((log_change, deltas))
    np.seterr(divide="warn", invalid="warn")
    # only choose from taxa that are not zero in 95% of the samples
    n_samples = p.shape[0]
    s = (p < 1e-4).sum(axis=0, keepdims = True) 
    taxa_list = np.where(s < 0.95*n_samples)
    # pick taxon with smallest change in log proportion
    min_idx = -1
    min_var = np.inf
    # ntaxa = log_change.shape[1]
    # for i in range(ntaxa):
    for i in taxa_list[1]:
        if not np.all(np.isfinite(log_change[:,i])):
            continue
        var = np.var(log_change[:,i])
        if var < min_var:
            min_idx = i
            min_var = var

    if min_idx == -1:
        print("Error: no valid denominator found", file=sys.stderr)
        exit(1)

    return min_idx


def construct_alr(P, denom, pseudo_count=1e-3):
    """Compute the additive log ratio transformation with a given
    choice of denominator. Assumes zeros have been replaced with
    nonzero values.
    """
    ALR = []
    ntaxa = P[0].shape[1]
    numer = np.array([i for i in range(ntaxa) if i != denom])
    for p in P:
        p = np.copy(p)
        p = (p + pseudo_count) / (p + pseudo_count).sum(axis=1,keepdims=True)
        p /= p.sum(axis=1, keepdims=True) 
        alr = (np.log(p[:,numer]).T - np.log(p[:,denom])).T
        ALR.append(alr)
    return ALR
