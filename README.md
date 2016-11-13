[![Build Status](https://travis-ci.org/casadi/python-optistack.png?branch=master)](https://travis-ci.org/casadi/python-optistack)

# optistack
The goal of this project is to provide a Yalmip-like wrapper around the Python interface of [CasADi](http://casadi.org), resulting in an easy and fast way to solve NLPs with [Ipopt](https://projects.coin-or.org/Ipopt) using exact sensitivities. 

Simple NLP:
```python
from optistack import *
x = optivar()
y = optivar()

nlp = optisolve((1-x)**2+100*(y-x**2)**2,[x**2+y**2<=1, x+y>=0])

print optival(x)
print optival(y)
```

Parametric NLP:
```python
from optistack import *
x=optivar(3,1)

# parameter, fixed during optimization
r=optipar()

r.setValue(1)
sol = optisolve(x[0]+x[1],[mtimes(x.T,x)==r])

print optival(x) # [-0.7071;-0.7071;0]

# adapt the value of the parameter
r.setValue(2)
sol.resolve()

print optival(x) # [-1;-1;0]
```

Installation:
 * Obtain CasADi for your version of Python:

Windows   |   Linux     |    Mac
----------|-------------|--------------
[Python 2.7](http://files.casadi.org/3.0.0/windows/casadi-py27-np1.9.1-v3.0.0.zip)  |    [Python 2.7](http://files.casadi.org/3.0.0/linux/casadi-py27-np1.9.1-v3.0.0.tar.gz) or later      | [Python 2.7](http://files.casadi.org/3.0.0/osx/casadi-py27-np1.9.1-v3.0.0.tar.gz) or later

 * Add the unzipped directory to the Python path (`import sys;sys.path.append(r"casadi_unzippeddir")`)
 * Obtain [Optistack](https://github.com/casadi/python-optistack/archive/master.zip) and unzip it.
 * Add the src directory to the Python path (`import sys;sys.path.append(r"optistack_unzippeddir")`)
 * Run the above example to verify that the install succeeded


Cheat sheet:

                          |  Yalmip                     | Optistack
------------------------- | --------------------------- | -----------------------------
Declare decision variable | `x = sdpvar(nrows,ncols)`   | `x = optivar(nrows,ncols)`
Set initial values        | `assign(x,x0)`              | `x.setInit(x0)`
Solve a problem           | `optimize([x==0],f)`        | `optisolve(f,[x==0])`
Obtain numeric results    | `value(x)`                  | `optival(x)`


More functionality:
`optivar` inherits from `casadi.MX` class. This means [all the usual operations](http://casadi.sourceforge.net/v3.0.0/api/html/d9/dc2/group__expression__tools.html) for `MX` are possible.


Note: there is also a [matlab version](http://optistack.casadi.org) available.



