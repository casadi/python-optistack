from optistack import *
from numpy import fabs

x = optivar()
y = optivar()

optisolve((x-1)**2+x)

print optival(x)

assert(max(fabs(optival(x)-0.5)<1e-9))

##

p = optipar()

p.setValue(100)

nlp = optisolve((1-x)**2+p*(y-x**2)**2,[x**2+y**2<=1, x+y>=0])

print optival(x)
print optival(y)

assert(max(fabs(optival(x)-0.7864151510041)<1e-9))
assert(max(fabs(optival(y)-0.617698307517294)<1e-9))
assert(max(fabs(optival(x**2+y**2)-1)<1e-7))

p.setValue(10)
nlp.resolve()

assert(max(fabs(optival(x)-0.788740490289645)<1e-9))
assert(max(fabs(optival(y)-0.614726302810607)<1e-9))
assert(max(fabs(optival(x**2+y**2)-1)<1e-7))

##
nlp = optisolve((1-x)**2+100*(y-x**2)**2,[x**2+y**2==1])

print optival(x)
print optival(y)

assert(max(fabs(optival(x)-0.7864151510041)<1e-8))
assert(max(fabs(optival(y)-0.617698307517294)<1e-8))
assert(max(fabs(optival(x**2+y**2)-1)<1e-7))

##
nlp = optisolve((1-x)**2+100*(y-x**2)**2,[0<=x, x<=0.5,0<=y, y<=0.5])

print optival(x)
print optival(y)

assert(max(fabs(optival(x)-0.5)<1e-8))
assert(max(fabs(optival(y)-0.25)<1e-8))

##
nlp = optisolve((1-x)**2+100*(y-x**2)**2,[1.5<=x, x<=2,1.5<=y, y<=2])

print optival(x)
print optival(y)

assert(max(fabs(optival(x)-1.5)<1e-8))
assert(max(fabs(optival(y)-2)<1e-8))
