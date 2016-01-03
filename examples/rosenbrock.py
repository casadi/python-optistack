from optistack import *
from numpy import fabs

x = optivar()
y = optivar()

nlp = optisolve((1-x)**2+100*(y-x**2)**2,[x**2+y**2<=1, x+y>=0])

print optival(x)
print optival(y)

assert(max(fabs(optival(x)-0.7864151510041)<1e-9))
assert(max(fabs(optival(y)-0.617698307517294)<1e-9))
assert(max(fabs(optival(x**2+y**2)-1)<1e-7))
