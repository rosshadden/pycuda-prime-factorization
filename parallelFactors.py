import sys

import nzmath.arith1 as arith1
import nzmath.factor.util as util
import nzmath.prime as prime

import pycuda.gpuarray as gpuarray
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule
import numpy


class ParallelTrialDivision (util.FactoringMethod):
    def __init__(self):
        util.FactoringMethod.__init__(self)

    def find(self, target, **options):
        return pTrialDivision(target)


kernel = SourceModule("""
    __global__ void factor(float *a){
        int idx = threadIdx.x + threadIdx.y*4;
        a[idx] *= 2;
    }
""")


def pTrialDivision(n, **options):
    if n < 1000000:
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


def main(n, **options):
    parallelTrialDivision = ParallelTrialDivision()
    return parallelTrialDivision.factor(n, **options)


try:
    n = int(sys.argv[1:][0])
except:
    n = 100

# test
print n
print main(n)
