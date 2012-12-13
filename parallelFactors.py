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
    __global__ void factor(long long n, int *value, int *primes){
        int i = threadIdx.x;
        int p = primes[i];

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
    numThreads = 384
    function = kernel.get_function('factor')

    values = numpy.zeros(len(allPrimes), numpy.int32)
    primes = numpy.copy(allPrimes).astype(numpy.int32)
    values_gpu = cuda.mem_alloc(values.nbytes)
    primes_gpu = cuda.mem_alloc(primes.nbytes)

    while True:
        result = numpy.array([], numpy.int32)
        # for t in range(0, int(numTimes)):
        cuda.memcpy_htod(values_gpu, values)
        cuda.memcpy_htod(primes_gpu, primes)

        function(numpy.int64(n), values_gpu, primes_gpu, block=(numThreads, 1, 1))
        currentResult = numpy.empty_like(values)
        cuda.memcpy_dtoh(currentResult, values_gpu)
        result = numpy.append(result, currentResult)

        factor = min2(result, 1)
        if factor == None:
            break
        factors.append(factor)
        n /= factor

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

n0 = n

result = factor(n)
print n, 'factored in', method
print result

total = 1
for r in result:
    total *= r
print 'Total multiplication of factors (for verification):\t', total
print '\tinput == output:\t', n0 == total
