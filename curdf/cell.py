import torch
from torch import Tensor


def cell_tensor(cell: Tensor | list | tuple, device, dtype) -> Tensor:
    """
    Convert an input cell into a torch tensor shaped ``(1, 3, 3)``.

    Parameters
    ----------
    cell
        Cell matrix shaped ``(3, 3)`` (triclinic allowed).
    device
        Torch device for the returned tensor.
    dtype
        Torch dtype for the returned tensor.

    Returns
    -------
    Tensor
        Cell tensor with leading batch dimension ``(1, 3, 3)``.
    """
    cell_t = torch.as_tensor(cell, device=device, dtype=dtype)
    if cell_t.shape != (3, 3):
        raise ValueError(f"cell must be (3,3); got {tuple(cell_t.shape)}")
    return cell_t.unsqueeze(0)


def pbc_tensor(pbc: Tensor | list | tuple, device) -> Tensor:
    """
    Convert periodic boundary flags to a boolean tensor shaped ``(1, 3)``.

    Parameters
    ----------
    pbc
        Iterable or tensor with three boolean flags.
    device
        Torch device for the returned tensor.

    Returns
    -------
    Tensor
        Boolean PBC tensor shaped ``(1, 3)``.
    """
    pbc_t = torch.as_tensor(pbc, device=device, dtype=torch.bool)
    if pbc_t.shape != (3,):
        raise ValueError(f"pbc must be (3,); got {tuple(pbc_t.shape)}")
    return pbc_t.unsqueeze(0)


def cell_volume(cell: Tensor) -> float:
    """
    Compute the volume of a cell tensor shaped ``(1, 3, 3)``.

    Parameters
    ----------
    cell
        Cell tensor with leading batch dimension ``(1, 3, 3)``.

    Returns
    -------
    float
        Absolute value of the cell determinant (volume).
    """
    if cell.shape != (1, 3, 3):
        raise ValueError(f"cell must be (1,3,3); got {tuple(cell.shape)}")
    return torch.det(cell[0]).abs().item()
