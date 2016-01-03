from OptimizationObject import OptimizationObject
import casadi as C

def optival( *args ):
    if len(args)==1 and isinstance(args[0],OptimizationObject):
      return args[0].optival().toArray()
    symbols = OptimizationObject.get_primitives(args)
    
    hassymbols = False;
    symbolsx = [C.MX.sym('dummy')]; # bug in casadi typemaps: [] does not work
    if 'x' in symbols:
       symbolsx = symbols["x"]
       hassymbols = True
    
    f = C.MXFunction('f',symbolsx,args);
    
    if hassymbols:
        for i,e in enumerate(symbolsx):
           f.setInput(optival(e),i)
    
    f.evaluate()
    
    if len(args)==1:
      return f.getOutput().toArray()
    else:
      return [f.getOutput(i).toArray() for i in range(f.nOut())]

