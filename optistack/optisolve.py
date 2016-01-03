from OptimizationObject import OptimizationObject
from sort_constraints import sort_constraints
import casadi as C
import numpy as np

class optisolve:
    
    def  __init__(self,objective,*args):
        """
               optisolve(objective)
               optisolve(objective,constraints)
        """
        if len(args)>=1:
            constraints = args[0]
        else:
            constraints = []
        options = dict()
        if len(args)>=2:
            options = args[1]
        
        
        if not isinstance(constraints,list):
            raise Exception("Constraints must be given as list: [x>=0,y<=0]")
            

        [ gl_pure, gl_equality] = sort_constraints( constraints )
        symbols = OptimizationObject.get_primitives([objective]+gl_pure)

        # helper functions for 'x'
        X = C.veccat(symbols["x"])
        helper = C.MXFunction('helper',[X],symbols["x"])

        helper_inv = C.MXFunction('helper_inv',symbols["x"],[X])

        # helper functions for 'p' if applicable
        if 'p' in symbols:
          P = C.veccat(symbols["p"])

          self.Phelper_inv = C.MXFunction('Phelper_inv',symbols["p"],[P])
          
        else:
          P = C.MX.sym('p',0,1)

        if len(gl_pure)>0:
            g_helpers = [];
            for p in gl_pure:
               g_helpers.append(C.MX.sym('g',p.sparsity())) 
            
            G_helpers = C.veccat(g_helpers)

            self.Ghelper = C.MXFunction('Ghelper',[G_helpers],g_helpers)

            self.Ghelper_inv = C.MXFunction('Ghelper_inv',g_helpers,[G_helpers])
        
        codegen = False;
        if 'codegen' in options:
            codegen = options["codegen"]
            del options["codegen"]
        
        opt = dict()
        if codegen:
            opt["jit"] = True;
            opt["jit_options"] = {'flags':['-O3']}
        
        gl_pure_v = C.MX()
        if len(gl_pure)>0:
           gl_pure_v = C.veccat(gl_pure)

        if objective.isvector() and objective.numel()>1:
            F = C.vec(objective)
            objective = 0.5*C.inner_prod(F,F)
            FF = C.MXFunction('nlp',[X,P], [F])
            JF = FF.jacobian()
            J_out = JF([X,P])
            J = J_out[0].T;
            H = C.mul(J,J.T)
            sigma = C.MX.sym('sigma')
            Hf = C.MXFunction('H',C.hessLagIn(x=X,p=P,lam_f=sigma),C.hessLagOut(hess=sigma*H),opt)
            if "expand" in options and options["expand"]:
               Hf = C.SXFunction(Hf)
            options["hess_lag"] = Hf
        
        nlp = C.MXFunction('nlp',C.nlpIn(x=X,p=P), C.nlpOut(f=objective,g=gl_pure_v),opt)

        self.solver = C.NlpSolver('solver','ipopt', nlp, options)

        # Save to class properties
        self.symbols      = symbols
        self.helper       = helper
        self.helper_inv   = helper_inv
        self.gl_equality  = gl_equality
        
        self.resolve()
          
    def jacspy(self):
        sparse(DMatrix(self.solver.jacG().output(),1)).spy()
    
    def resolve(self):
        # recall from class properties
        symbols      = self.symbols
        helper       = self.helper
        helper_inv   = self.helper_inv
        gl_equality  = self.gl_equality
        if len(gl_equality)>0:
            # compose lbg
            for i,e in enumerate(gl_equality):
                if e:
                    self.Ghelper_inv.setInput(0,i)
                else:
                    self.Ghelper_inv.setInput(-np.inf,i)
            self.Ghelper_inv.evaluate()
            self.solver.setInput(self.Ghelper_inv.getOutput(),'lbg')
            self.solver.setInput(0,'ubg')

        # compose lbx
        for i,e in enumerate(symbols["x"]):
          helper_inv.setInput(e.lb,i)

        helper_inv.evaluate()
        self.solver.setInput(helper_inv.getOutput(),'lbx')

        # compose x0
        for i,e in enumerate(symbols["x"]):
          helper_inv.setInput(e.init,i)

        helper_inv.evaluate();
        self.solver.setInput(helper_inv.getOutput(),'x0')

        # compose ubx
        for i,e in enumerate(symbols["x"]):
          helper_inv.setInput(e.ub,i)

        helper_inv.evaluate()
        self.solver.setInput(helper_inv.getOutput(),'ubx')

        if 'p' in symbols:
            # compose p0
            for i,e in enumerate(symbols["p"]):
              self.Phelper_inv.setInput(e.value,i)

            self.Phelper_inv.evaluate()
            self.solver.setInput(self.Phelper_inv.getOutput(),'p')

        self.solver.evaluate()

        helper.setInput(self.solver.getOutput('x'))
        helper.evaluate()

        for i,e in enumerate(symbols["x"]):
          e.setValue(helper.getOutput(i))
