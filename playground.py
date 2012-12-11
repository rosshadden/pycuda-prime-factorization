from __future__ import generators
import sys, math, fractions
import nzmath.arith1 as arith1
import nzmath.factor.util as util
import nzmath.gcd as gcd
import nzmath.prime as prime

## {{{ http://code.activestate.com/recipes/117119/ (r2)
# Sieve of Eratosthenes
# David Eppstein, UC Irvine, 28 Feb 2002


def eratosthenes():
    '''Yields the sequence of prime numbers via the Sieve of Eratosthenes.'''
    D = {}  # map composite integers to primes witnessing their compositeness
    q = 2   # first integer to test for primality
    while 1:
        if q not in D:
            yield q        # not marked composite, must be prime
            D[q * q] = [q]   # first multiple of q not already marked
        else:
            for p in D[q]:  # move each witness to its next multiple
                D.setdefault(p + q, []).append(p)
            del D[q]       # no longer need D[q], free memory
        q += 1
## end of http://code.activestate.com/recipes/117119/ }}}


def factor(n):
    """Return a list of the prime factors for a natural number."""
    if n == 1:
        return [1]
    primes = sieve(int(n ** 0.5) + 1)
    prime_factors = []

    for p in primes:
        if p * p > n:
            break
        while n % p == 0:
            prime_factors.append(p)
            n //= p
    if n > 1:
        prime_factors.append(n)

    return prime_factors


class ParallelTrialDivision (util.FactoringMethod):
    def __init__(self):
        util.FactoringMethod.__init__(self)

    def find(self, target, **options):
        return pTrialDivision(target)


def parallelTrialDivision(n, **options):
    options['return_type'] = 'list'
    options['need_sort'] = True
    parallelTrialDivision = ParallelTrialDivision()
    return parallelTrialDivision.factor(n, **options)


def pTrialDivision(n, **options):
    # if 'start' in options and 'stop' in options:
        # if 'step' in options:
            # trials = range(options['start'], options['stop'], options['step'])
        # else:
            # trials = range(options['start'], options['stop'])
    if 'iterator' in options:
        trials = options['iterator']
    elif n < 1000000:
        trials = prime.generator_eratosthenes(arith1.floorsqrt(n))
    else:
        trials = prime.generator()

    limit = arith1.floorsqrt(n)
    for p in trials:
        if limit < p:
            break
        if 0 == n % p:
            return p
    return 1


def parallelPMinusOne(n):
    # initialize
    x = y = 2
    primes = []
    B = 10000

    for q in prime.generator():
        primes.append(q)
        if q > B:
            if gcd.gcd(x - 1, n) == 1:
                return 1
            x = y
            break
        q1 = q
        l = B // q
        while q1 <= l:
            q1 *= q
        x = pow(x, q1, n)
        if len(primes) >= 20:
            if gcd.gcd(x - 1, n) == 1:
                primes, y = [], x
            else:
                x = y
                break

    for q in primes:
        q1 = q
        while q1 <= B:
            x = pow(x, q, n)
            g = gcd.gcd(x - 1, n)
            if g != 1:
                if g == n:
                    return 1
                return g
            q1 *= q


def factors(n):
    result = []
    for i in range(2, n + 1):  # test all integers between 2 and n
        s = 0
        while n / float(i) == math.floor(n / float(i)):  # is n/i an integer?
            n /= float(i)
            s += 1
        if s > 0:
            for k in range(s):
                result.append(i)  # i is a prime factor s times
            if n == 1:
                return result
    return result


def factorize(n):
    "Returns all the prime factors of a positive integer"
    factors = []
    d = 2
    while (n > 1):
        while n % d is 0:
            factors.append(d)
            n /= d
        d = d + 1
        if d * d > n:
            if n > 1:
                factors.append(n)
            break
    return factors


def sieve(n):
    '''Get all primes up to n.'''
    n = int(n)
    if n < 2:
        return []
    sieve = range(n)
    sieve[1] = 0
    root = n ** 0.5
    index = 0
    while index <= root:
        if sieve[index]:
            i = index ** 2
            while i < n:
                sieve[i] = 0
                i += index
        index += 1
    return [x for x in sieve if x]

try:
    n = int(sys.argv[1:][0])
except:
    n = 100

# test
print n
print factor(n)
