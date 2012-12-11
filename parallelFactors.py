import sys

import pycuda.gpuarray as gpuarray
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
    __global__ void factor(int n, int *value){
        int i = threadIdx.x;
        int p = value[i];

        if(n % p == 0){
            value[i] = p;
            //value[i] = n / p;
        }else{
            value[i] = 9999;
        }
    }
""")


def factor(n):
    """Return a list of the prime factors for a natural number."""
    if n == 1:
        return [1]
    primes = quadradticSieve(int(n ** 0.5) + 1)
    factors = []

    a = numpy.copy(primes)
    a = a.astype(numpy.int32)
    a_gpu = cuda.mem_alloc(a.nbytes)
    cuda.memcpy_htod(a_gpu, a)
    function = kernel.get_function('factor')

    # for p in primes:
        # if p * p > n:
        #     break
        # while n % p == 0:
        #     factors.append(p)
        #     n //= p
    while True:
        function(numpy.int32(n), a_gpu, block=(128, 1, 1))
        a_copy = numpy.empty_like(a)
        cuda.memcpy_dtoh(a_copy, a_gpu)
        factor = min(a_copy)
        if factor == 9999 or factor == 1 or factor == 0:
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


try:
    n = int(sys.argv[1:][0])
except:
    n = 100

# test
print n
print factor(n)
