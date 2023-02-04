import torch

def evp(mor: torch.Tensor, basis: torch.Tensor) -> torch.Tensor:
    _mor = mor.clone().detach()
    _p = (_mor**basis).prod(dim=1).sum()
    return _p.clone().detach()

def evdpdmor(mor: torch.Tensor, basis: torch.Tensor) -> torch.Tensor:
    _mor = mor.clone().detach()
    _mor.requires_grad = True
    _p = (_mor**basis).prod(dim=1).sum()
    _p.backward(torch.ones_like(_p))
    grad = _mor._grad.clone().detach()
    return grad
