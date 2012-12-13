import sys
import math

import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy


def quadradticSieve(n):
    """Get all primes up to n."""
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


kernel = SourceModule("""
    __global__ void factor(long long n, int *value){
        int i = threadIdx.x;
        int p = value[i];

        if(n % p == 0){
            value[i] = p;
        }else{
            value[i] = 0;
        }
    }
""")


def factorParallel(n):
    """Return a list of the prime factors for a natural number."""

    def min2(list, bound=0):
        for item in list:
            if item > bound:
                return item
        return None

    if n == 1:
        return [1]
    factors = []

    allPrimes = quadradticSieve(int(n ** 0.5) + 1)
    function = kernel.get_function('factor')

    logged = False
    while True:
        limit = min(len(allPrimes), 384)
        numTimes = math.ceil(len(allPrimes) / (1.0 * limit))

        result = numpy.array([], numpy.int32)
        for t in range(0, int(numTimes)):
            primes = numpy.copy(allPrimes[t * limit:min(t * limit + limit, len(allPrimes))]).astype(numpy.int32)
            primes_gpu = cuda.mem_alloc(primes.nbytes)
            cuda.memcpy_htod(primes_gpu, primes)

            function(numpy.int64(n), primes_gpu, block=(limit, 1, 1))
            currentResult = numpy.empty_like(primes)
            cuda.memcpy_dtoh(currentResult, primes_gpu)
            result = numpy.append(result, currentResult)
            if min2(result, 1) != None:
                break

        factor = min2(result, 1)
        if factor == None:
            break
        factors.append(factor)
        n /= factor

        if not logged:
            logged = True
            print 'Using', len(result), 'threads in parallel.'
        print min2(result, 1), '\t', n, '\t', numpy.int64(n)

    if n > 1:
        factors.append(n)

    return factors


def factorSerial(n):
    """Return a list of the prime factors for a natural number."""
    if n == 1:
        return [1]
    primes = quadradticSieve(int(n ** 0.5) + 1)
    factors = []

    for p in primes:
        if p * p > n:
            break
        while n % p == 0:
            factors.append(p)
            n //= p

    if n > 1:
        factors.append(n)

    return factors


if len(sys.argv) == 3:
    method = sys.argv[1]
    n = int(sys.argv[2])
elif len(sys.argv) == 2:
    n = int(sys.argv[1])
    method = 'parallel'
else:
    n = 100
    method = 'parallel'

if method == 'serial':
    factor = factorSerial
else:
    factor = factorParallel

# test
result = factor(n)
print n, 'factored in', method
print result

total = 1
for r in result:
    total *= r
print 'Total multiplication of factors (for verification):\t', total


# use stuff around 2147483647.
