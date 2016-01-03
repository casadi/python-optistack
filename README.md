[![Build Status](https://travis-ci.org/casadi/python-optistack.png?branch=master)](https://travis-ci.org/casadi/python-optistack)

# optistack
The goal of this project is to provide a Yalmip-like wrapper around the Python interface of [CasADi](http://casadi.org).

Example:
```mpython
x = optivar()
y = optivar()

nlp = optisolve((1-x)**2+100*(y-x**2)**2,[x**2+y**2<=1, x+y>=0])

print optival(x)
print optival(y)
```

Installation:
 * Obtain CasADi for your version of Python:

Windows   |   Linux     |    Mac
----------|-------------|--------------
[Python 2.7](http://files.casadi.org/2.4.2/windows/casadi-py27-np1.9.1-v2.4.2.zip)  |    [Python 2.7](http://files.casadi.org/2.4.2/linux/casadi-py27-np1.9.1-v2.4.2.tar.gz) or later      | [Python 2.7](http://files.casadi.org/2.4.2/osx/casadi-py27-np1.9.1-v2.4.2.tar.gz) or later

 * Add the unzipped directory to the Python path (`import sys;sys.path.append('casadi_unzippeddir')`)
 * Obtain [Optistack](https://github.com/casadi/python-optistack/archive/master.zip) and unzip it.
 * Add the src directory to the Python path (`import sys;sys.path.append('optistack_unzippeddir')`)
 * Run the above example to verify that the install succeeded


Cheat sheet:

                          |  Yalmip                     | Optistack
------------------------- | --------------------------- | -----------------------------
Declare decision variable | `x = sdpvar(nrows,ncols)`   | `x = optivar(nrows,ncols)`
Set initial values        | `assign(x,x0)`              | `x.setInit(x0)`
Solve a problem           | `optimize([x==0],f)`        | `optisolve(f,[x==0])`
Obtain numeric results    | `value(x)`                  | `optival(x)`


More functionality:
`optivar` inherits from `casadi.MX` class. This means [all the usual operations](http://casadi.sourceforge.net/v2.4.1/api/html/d9/dc2/group__expression__tools.html) for `MX` are possible.






