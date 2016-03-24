from OptimizationObject import OptimizationObject
import casadi as C
import numpy as np

def optival( *args ):
    if len(args)==1 and isinstance(args[0],OptimizationObject):
      return np.array(args[0].optival())
    symbols = OptimizationObject.get_primitives(args)
    
    hassymbols = False;
    symbolsx = [C.MX.sym('dummy')]; # bug in casadi typemaps: [] does not work
    if 'x' in symbols:
       symbolsx = symbols["x"]
       hassymbols = True
    
    f = C.Function('f',symbolsx,args);
    f_inputs = [];
    
    if hassymbols:
        for i,e in enumerate(symbolsx):
           f_inputs.append(optival(e))
    
    out = f.call(f_inputs)
    
    if len(args)==1:
      return np.array(out[0])
    else:
      return [np.array(out[i]) for i in range(f.nOut())]

