import torch
from torch import Tensor


def cell_tensor(cell: Tensor | list | tuple, device, dtype) -> Tensor:
    """
    Convert an input (3,3) cell array into a torch tensor shaped (1,3,3).
    Accepts triclinic cells. Does not assume orthorhombic.
    """
    cell_t = torch.as_tensor(cell, device=device, dtype=dtype)
    if cell_t.shape != (3, 3):
        raise ValueError(f"cell must be (3,3); got {tuple(cell_t.shape)}")
    return cell_t.unsqueeze(0)


def pbc_tensor(pbc: Tensor | list | tuple, device) -> Tensor:
    """
    Convert PBC flags to shape (1,3) boolean tensor.
    """
    pbc_t = torch.as_tensor(pbc, device=device, dtype=torch.bool)
    if pbc_t.shape != (3,):
        raise ValueError(f"pbc must be (3,); got {tuple(pbc_t.shape)}")
    return pbc_t.unsqueeze(0)


def cell_volume(cell: Tensor) -> float:
    """
    Compute volume from a (1,3,3) cell tensor (triclinic supported).
    """
    if cell.shape != (1, 3, 3):
        raise ValueError(f"cell must be (1,3,3); got {tuple(cell.shape)}")
    return torch.det(cell[0]).abs().item()
