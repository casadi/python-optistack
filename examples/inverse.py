import numpy as np
from optistack import *
import casadi as C

A = optivar(3,3)

# inverse
B = np.array([[1, 0, -2],[3, 1, -1],[2, 1, 0]])

optisolve(C.sum_square(mul(A,B)-np.eye(3)))

As = optival(A)

assert(np.max(np.max(np.fabs(As-np.linalg.inv(B)))<1e-7))

A = optivar(4,3);
# moore-penrose pseudo-inverse
B = np.array([[1, 0, -2, 3],[ 3, 1, -1, 1],[2, 1, 0, 9]])

optisolve(C.sum_square(mul(A,B)-np.eye(4)))

As = optival(A)

assert(np.max(np.max(np.fabs(As-np.linalg.pinv(B)))<1e-7))






