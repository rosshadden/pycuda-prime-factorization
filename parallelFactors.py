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

    factor = kernel.get_function('factor')

    if n == 1:
        return [1]
    factors = []

    sqrtN = int(n ** 0.5) + 1
    allPrimes = quadradticSieve(sqrtN)
    numPrimes = len(allPrimes)
    numThreads = 384

    values = numpy.zeros(numPrimes, numpy.int32)
    primes = numpy.copy(allPrimes).astype(numpy.int32)
    values_gpu = cuda.mem_alloc(values.nbytes)
    primes_gpu = cuda.mem_alloc(primes.nbytes)
    cuda.memcpy_htod(values_gpu, values)
    cuda.memcpy_htod(primes_gpu, primes)

    result = numpy.zeros(numPrimes, numpy.int32)

    while True:
        factor(numpy.int64(n), values_gpu, primes_gpu, block=(numThreads, 1, 1))
        cuda.memcpy_dtoh(result, values_gpu)

        prime = min2(result, 1)
        if prime == None:
            break
        factors.append(prime)
        n /= prime

        sqrtN = int(n ** 0.5) + 1
        numPrimesUnderN = math.ceil(sqrtN / math.log(sqrtN))
        primes = primes[0:numPrimesUnderN]
        cuda.memcpy_htod(primes_gpu, primes)

    print len(allPrimes)

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
