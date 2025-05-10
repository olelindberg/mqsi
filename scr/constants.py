from dataclasses import dataclass

@dataclass(frozen=True)
class constants:
    NODE_DOFS: int = 6