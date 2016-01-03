from optistack import *
from casadi import *

# In this example, we fit a nonlinear model to measurements
#
# This example uses more advanced constructs than the vdp* examples:
# Since the number of control intervals is potentially very large here,
# we use memory-efficient Map and MapAccum, in combination with
# codegeneration.
#
# We will be working with a 2-norm objective:
# || y_measured - y_simulated ||_2^2
#
# This form is well-suited for the Gauss-Newton Hessian approximation.

########### SETTINGS #####################
N = 10000  # Number of samples
fs = 610.1 # Sampling frequency [hz]

param_truth = np.array([5.625e-6,2.3e-4,1,4.69])
param_guess = np.array([5,2,1,5])
scale = np.array([1e-6,1e-4,1,1])

############ MODELING #####################
y  = MX.sym('y')
dy = MX.sym('dy')
u  = MX.sym('u')

states = vertcat([y,dy])
controls = u

M = optivar()
c = optivar()
k = optivar()
k_NL = optivar()

params = vertcat([M,c,k,k_NL])

rhs = vertcat([dy, (u-k_NL*y**3-k*y-c*dy)/M])

# Form an ode function
ode = MXFunction('ode',[states,controls,params],[rhs])

############ Creating a simulator ##########
N_steps_per_sample = 10
dt = 1/fs/N_steps_per_sample

# Build an integrator for this system: Runge Kutta 4 integrator
k1 = ode([states,controls,params])
k2 = ode([states+dt/2.0*k1[0],controls,params])
k3 = ode([states+dt/2.0*k2[0],controls,params])
k4 = ode([states+dt*k3[0],controls,params])

states_final = states+dt/6.0*(k1[0]+2*k2[0]+2*k3[0]+k4[0])

# Create a function that simulates one step propagation in a sample
one_step = MXFunction('one_step',[states, controls, params],[states_final])

X = states
for i in range(N_steps_per_sample):
    Xn = one_step([X, controls, params])
    X = Xn[0]

# Create a function that simulates all step propagation on a sample
one_sample = MXFunction('one_sample',[states, controls, params], [X])

# speedup trick: expand into scalar operations
one_sample = one_sample.expand()

############ Simulating the system ##########

all_samples = one_sample.mapaccum('all_samples', N)

# Choose an excitation signal
u_data = 0.1*np.random.rand(N,1)

x0 = DMatrix([0,0])
all_samples_out = all_samples([x0, u_data, repmat(param_truth,1,N) ])
X_measured = all_samples_out[0]

y_data = X_measured[0,:].T

############ Identifying the simulated system: single shooting strategy ##########

# Note, it is in general a good idea to scale your decision variables such
# that they are in the order of ~0.1..100
all_samples_out = all_samples([x0, u_data, repmat(params*scale,1,N) ])
X_symbolic = all_samples_out[0]

e = y_data-X_symbolic[0,:].T

M.setInit(param_guess[0])
c.setInit(param_guess[1])
k.setInit(param_guess[2])
k_NL.setInit(param_guess[3])

options = {}
options["codegen"] = True

print 'Single shooting...'

# Hand in a vector objective -> interpreted as 2-norm
# such t hat Gauss-Newton can be performed
optisolve(e,[],options)
optival(M)*1e-6
optival(c)*1e-4
optival(k)
optival(k_NL)

print sqrt(np.mean(optival(e)**2))
assert(sqrt(np.mean(optival(e)**2))<1e-10)

############ Identifying the simulated system: multiple shooting strategy ##########

print 'Multiple shooting...'

X = optivar(2, N)

params_scale = vertcat([1e-6*M,c*1e-4,k,k_NL])

Xn = one_sample.map([X, u_data.T, params_scale])
Xn = Xn[0]

# gap-closing constraints
gaps = Xn[:,:-1]-X[:,1:]
g = gaps == 0

e = (y_data-Xn[0,:].T).T

M.setInit(5)
c.setInit(2.3)
k.setInit(1)
k_NL.setInit(4)

x0 = horzcat([y_data, vertcat([np.diff(y_data.T).T*fs,0])])
X.setInit(x0.T)

options = {}
options["codegen"] = True

# Hand in a vector objective -> interpreted as 2-norm 
# such that Gauss-Newton can be performed
optisolve(e,[g],options);

optival(M)*1e-6
optival(c)*1e-4
optival(k)
optival(k_NL)

# rms(optival(e))
print sqrt(np.mean(optival(e)**2)) # rms is part of the signal toolbox
assert(sqrt(np.mean(optival(e)**2))<1e-10)
