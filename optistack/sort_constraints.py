import casadi as C

def sort_constraints( gl ):
    """
     Rewrites and determines nature of constraints, either g(x)<=0 or g(x)==0.
     A user may write x>=y where x and y are variables.
     In the `gl_pure` output, everything is brought to the left hand side
     Parameters
     ----------
     gl : list of constraints, optional
     Returns
     -------
     gl_pure : list of constraints in standard form
     The constraints are rewritten as g(x)<=0 or g(x)==0
     gl_equality : list of bools
     For each entry in `gl_pure`, this list contains a boolean.
    """
    
    gl_pure = []
    gl_equality = []
    for g in gl:
        if g.is_op(C.OP_LE) or g.is_op(C.OP_LT):
          gl_pure.append(g.dep(0) - g.dep(1))
          gl_equality.append(False)
        elif g.is_op(C.OP_EQ):
          gl_pure.append(g.dep(0) - g.dep(1))
          gl_equality.append(True)
        else:
            raise Exception('Constraint type unkown. Use ==, >= or <= .');
    return [ gl_pure, gl_equality]
