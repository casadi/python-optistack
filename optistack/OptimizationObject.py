import casadi as C
import numpy as np

class OptimizationObject(C.MX):
    mapping = {}

    def __init__(self,shape,name):
       if isinstance(shape,int):
           shape = [shape, 1]
       C.MX.__init__(self,C.MX.sym(name,C.Sparsity.dense(*shape)))
       mymapping = OptimizationObject.mapping
       mymapping[hash(self)] = self
       self.value = np.nan
    
    def optival(self):
        return self.value

    def setValue(self,value):
        self.value = value
    
    def  reshape(self,a,b):
       if self.size1() == a and self.size2()==b:
          return self
       else:
          return casadi.reshape(self,(a,b))

    @staticmethod
    def get_primitives( el, *args):
            """
              Out of a list of expression, retrieves all primitive expressions
              The result is sorted into a dictionary with the key originating
              from the 'shorthand' class attribute of OptimzationObject subclasses
            """
            if len(args)==0:
               dep = True
            else:
               dep = args[0]
            
            try:
                vars = C.symvar(C.veccat(el))
            except:
                vars = {}
            
            syms = dict()

            for i,v in enumerate(vars):
                mymapping = OptimizationObject.mapping
                vobj = mymapping[hash(v)]
                symkey = vobj.shorthand
                if symkey in syms:
                    syms[symkey] = syms[symkey] + [vobj]
                else:
                    syms[symkey] = [vobj]

            return syms
