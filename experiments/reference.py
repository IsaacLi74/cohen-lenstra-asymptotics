"""Proved negative-discriminant torsion coefficient: BST, Theorem 4."""
from functools import lru_cache
import math

import numpy as np
from scipy.special import gamma, zeta

PRIME_BOUND = 1_000_000
SOURCE = "https://arxiv.org/pdf/1005.0672"


@lru_cache(maxsize=None)
def torsion_coefficient(prime_bound=PRIME_BOUND):
    """Euler product truncated at prime_bound; this is not a fitted value."""
    if prime_bound < 2:
        raise ValueError("Prime bound must be at least 2")
    sieve = np.ones(prime_bound+1, dtype=bool)
    sieve[:2] = False
    for p in range(2, math.isqrt(prime_bound)+1):
        if sieve[p]:
            sieve[p*p::p] = False
    primes = np.flatnonzero(sieve).astype(float)
    product = np.prod(1-(np.cbrt(primes)+1)/(primes*(primes+1)))
    return float(8*math.sqrt(3)*zeta(1/3)/(5*gamma(2/3)**3)*product)
