from OptimizationObject import OptimizationObject

class optipar(OptimizationObject):
    """
     Create an optimization parameter,
     something that is not optimised for.
    
      x = optipar()                   Scalar
      x = optipar(n)                  Column vector of length n
      x = optipar(n,m)                Matrix of shape nxn
      x = optipar(n,m,name)           A name for printing the variable
    """

    shorthand = 'p'
    
    def __init__(self,*args):
        shape = 1
        if len(args)==1:
           shape = args[0]
        elif len(args)>1:
           shape = [args[0],args[1]]

        if len(args)==3:
           name = args[2]
        else:
           name = 'p'
               
        OptimizationObject.__init__(self,shape,name)

