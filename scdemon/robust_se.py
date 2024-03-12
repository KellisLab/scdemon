import numpy as np
import scipy.sparse

from . import _core as core
from .progress_manager import ProgressManager, _interrupt_checker

def _robust_se_X(i:int, V:np.ndarray) -> np.ndarray:
    return core.py_robust_se_X(V[:, i].astype("f8"), V.astype("f8"))

def _robust_se(X:np.ndarray, V:np.ndarray, t_cutoff:float, abs_t:bool) -> scipy.sparse.csc_matrix:
    return core.py_robust_se(X.astype("f8"), V.astype('f8'), t_cutoff, abs_t)

def ols_resid(X:np.ndarray, Y:np.ndarray, beta:np.ndarray) -> np.ndarray:
    return core.py_ols_resid(X.astype("f8"), Y.astype("f8"), beta.astype("f8"))

def ols_beta(X:np.ndarray, Y:np.ndarray) -> np.ndarray:
    return core.py_ols_beta(X.astype("f8"), Y.astype("f8"))

def robust_prepare(U:np.ndarray, V:np.ndarray, B=None, n_components:int=None, min_norm:float=1e-5, return_U=False):
    if U.shape[1] != V.shape[0]:
        raise ValueError("Shapes must match")
    if n_components != None:
        if n_components <= 0:
            raise ValueError("N components must be > 0")
        U = U[:, :n_components]
        V = V[:n_components, :]
    if B is None:
        B = np.ones((U.shape[0], 1))
    if U.shape[0] != B.shape[0]:
        raise ValueError("U and B are not the correct shape")
    V1 = np.nan * np.ones_like(V)
    bad_flag = np.linalg.norm(V, axis=0, ord=2) < min_norm
    V_svd = np.linalg.svd(V[:, ~bad_flag], full_matrices=False)
    U = U @ V_svd.U
    qr = np.linalg.qr(ols_resid(X=B, Y=U, beta=ols_beta(X=B, Y=U)))
    RS_svd = np.linalg.svd(qr.R @ np.diag(V_svd.S))
    V1[:, ~bad_flag] = np.diag(RS_svd.S) @ RS_svd.Vh @ V_svd.Vh
    dof = B.shape[0] - B.shape[1]
    if return_U:
        return V1, dof, qr.Q @ RS_svd.U
    else:
        return V1, dof, None

def robust_se(U, V, B=None, t_cutoff:float=None, abs_t:bool=False, nominal_p_cutoff:float=0.05) -> scipy.sparse.csc_matrix:
    """
    U: U from SVD
    V: V\Sigma from SVD
    B: Batch matrix from pd.get_dummies, or just intercept.
    """
    if U.shape[1] != V.shape[0]:
        raise ValueError("U and V must have compatible dimensions : %d != %d" % (U.shape[1], V.shape[0]))
    if B is None:
        B = np.ones((U.shape[0], 1), dtype="f8")
    elif U.shape[0] != B.shape[0]:
        raise ValueError("U and B must have compatible dimensions : %d != %d" % (U.shape[0], B.shape[0]))
    UpB = ext.py_ols_beta(U.astype("f8"), B.astype("f8"))
    UpU = ext.py_ols_beta(U.astype("f8"), U.astype("f8"))
    if t_cutoff is None:
        import scipy.stats
        t_cutoff = scipy.stats.t.isf(nominal_p_cutoff * V.shape[1]**-2,
                                     B.shape[0] - B.shape[1])
    return _robust_se(V, V, t_cutoff=t_cutoff, abs_t=abs_t)


