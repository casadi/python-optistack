import numpy as np
from optistack import *
import time

n = 10
m = 10000

A = np.random.rand(m,n)
b = np.random.rand(m,1)

x = optivar(n,1,'x')

# This example shows a least-squares computation:
# min_x  || F(x) ||^2_2
#
# with F(x) = Ax+b
#

# A naive way to do this is:
t0 = time.time()
e = mul(A,x)+b
optisolve(0.5*mul(e.T,e))
print time.time()-t0

x1 = optival(x)

# Solving with Gauss-Newton: pass vector as objective
t0 = time.time()
e = mul(A,x)+b
optisolve(e)
print time.time()-t0

x2 = optival(x)

# Solutions are the same
assert(max(np.fabs(x1-x2))<1e-10)


