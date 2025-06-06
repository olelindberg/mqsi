from dataclasses import dataclass

@dataclass(frozen=True)
class Constants:
    NODE_DOFS: int = 6

def assign_constraints(x0, bc_dof, bc_value):
    for i in range(len(bc_dof)):
        for j in range(len(bc_dof[i])):
            k = i*Constants.NODE_DOFS + bc_dof[i][j]
            x0[k] = bc_value[i][j]
    return x0

def assign_constraints_grad(x0, bc_dof):
    for i in range(len(bc_dof)):
        for j in range(len(bc_dof[i])):
            k = i*Constants.NODE_DOFS + bc_dof[i][j]
            x0[k] = 0.0
    return x0