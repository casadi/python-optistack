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
        X = C.veccat(*symbols["x"])
        helper = C.Function('helper',[X],symbols["x"])

        helper_inv = C.Function('helper_inv',symbols["x"],[X])

        # helper functions for 'p' if applicable
        if 'p' in symbols:
          P = C.veccat(*symbols["p"])

          self.Phelper_inv = C.Function('Phelper_inv',symbols["p"],[P])
          
        else:
          P = C.MX.sym('p',0,1)

        if len(gl_pure)>0:
            g_helpers = [];
            for p in gl_pure:
               g_helpers.append(C.MX.sym('g',p.sparsity())) 
            
            G_helpers = C.veccat(*g_helpers)

            self.Ghelper = C.Function('Ghelper',[G_helpers],g_helpers)

            self.Ghelper_inv = C.Function('Ghelper_inv',g_helpers,[G_helpers])
        
        codegen = False;
        if 'codegen' in options:
            codegen = options["codegen"]
            del options["codegen"]
        
        opt = {}
        if codegen:
            options["jit"] = True
            opt["jit"] = True
        
        gl_pure_v = C.MX()
        if len(gl_pure)>0:
           gl_pure_v = C.veccat(*gl_pure)

        if objective.is_vector() and objective.numel()>1:
            F = C.vec(objective)
            objective = 0.5*C.dot(F,F)
            FF = C.Function('nlp',[X,P], [F])
            JF = FF.jacobian()
            J_out = JF.call([X,P])
            J = J_out[0].T;
            H = C.mtimes(J,J.T)
            sigma = C.MX.sym('sigma')
            Hf = C.Function('H',dict(x=X,p=P,lam_f=sigma,hess_gamma_x_x=sigma*C.triu(H)),['x', 'p', 'lam_f', 'lam_g'],['hess_gamma_x_x'],opt)
            if "expand" in options and options["expand"]:
               Hf = Hf.expand()
            options["hess_lag"] = Hf
        
        nlp = {"x":X,"p":P,"f":objective,"g":gl_pure_v}

        self.solver = C.nlpsol('solver','ipopt', nlp, options)

        # Save to class properties
        self.symbols      = symbols
        self.helper       = helper
        self.helper_inv   = helper_inv
        self.gl_equality  = gl_equality
        
        self.resolve()
          
    def jacspy(self):
        sparse(DM(self.solver.jacG().output(),1)).spy()
    
    def resolve(self):
        # recall from class properties
        symbols      = self.symbols
        helper       = self.helper
        helper_inv   = self.helper_inv
        gl_equality  = self.gl_equality
        solver_inputs = dict()
        if len(gl_equality)>0:
            Ghelper_inv_inputs = [[]]*self.Ghelper_inv.n_in()
            # compose lbg
            for i,e in enumerate(gl_equality):
                if e:
                    Ghelper_inv_inputs[i] = 0
                else:
                    Ghelper_inv_inputs[i] = -np.inf
            solver_inputs["lbg"] = self.Ghelper_inv.call(Ghelper_inv_inputs)[0]
            solver_inputs["ubg"] = 0

        helper_in_inputs = [[]]*helper_inv.n_in()
  
        # compose lbx
        solver_inputs["lbx"] = helper_inv.call([e.lb for e in symbols["x"]])[0]


        # compose x0
        solver_inputs["x0"] = helper_inv.call([e.init for e in symbols["x"]])[0]

        # compose ubx
        solver_inputs["ubx"] = helper_inv.call([e.ub for e in symbols["x"]])[0]

        if 'p' in symbols:
            # compose p0
            solver_inputs["p"] = self.Phelper_inv([e.value for e in symbols["p"]])[0]

        sol = self.solver.call(solver_inputs)

        out = helper.call([sol["x"]])

        for i,e in enumerate(symbols["x"]):
          e.setValue(out[i])
