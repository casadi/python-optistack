from numpy import *
from optistack import *
from rk4 import rk4

N = 100

# ode
ode = lambda x,u: vertcat([u - x[0], x[0]])

# state: speed
S = optivar(N+1,1,'S')
# state: position
P = optivar(N+1,1,'P')

# control: 
U = optivar(N,1,'U')

# a variable to denote the final time
tf = optivar(1,1,'tf')
tf.setInit(1)

# Construct list of all constraints
g = []

for k in range(N):
   xk      = vertcat([S[k],   P[k]  ])
   xk_plus = vertcat([S[k+1], P[k+1]])
   
   # shooting constraint
   xf = rk4(lambda x,u: tf*ode(x,u),1.0/N,xk,U[k])
   g.append(xk_plus==xf)


# path constraint
constr = lambda P: 1-sin(2*pi*P)/2
g.append(S <= constr(P))

# Initialize with speed 1.
S.setInit(1)

U.setLb(0)
U.setUb(1)

g+= [S[0]==0, P[0]==0, P[-1]==1]

optisolve(tf,g)

from pylab import *

plot(optival(S))
plot(optival(P))
plot(constr(optival(P)),'r--')
step(range(N),optival(U),'b')
show()
