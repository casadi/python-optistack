from OptimizationObject import OptimizationObject

class optivar(OptimizationObject):
    """
      Create an optimization variable
    
      x = optivar()                   Scalar
      x = optivar(n)                  Column vector of length n
      x = optivar(n,m)                Matrix of shape nxn
      x = optivar(n,m,name)           A name for printing the variable
    """
    shorthand = 'x'

    def setInit(self,v):
        self.init[:,:] = v
    def setLb(self,v):
        self.lb[:,:] = v
    def setUb(self,v):
        self.ub[:,:] = v
    def __init__(self,*args):
       shape = 1
       if len(args)==1:
           shape = args[0]
       elif len(args)>1:
           shape = [args[0],args[1]]
       
       if len(args)==3:
           name = args[2]
       else:
           name = 'x'
       
       import casadi
       OptimizationObject.__init__(self,shape,name)
       self.lb = -casadi.inf*casadi.DM.ones(self.sparsity())
       self.ub = casadi.inf*casadi.DM.ones(self.sparsity())
       self.init = casadi.DM.zeros(self.sparsity())
