import numpy as np
from scipy import sparse


#from .utils import ProgressManager, _interrupt_checker

def _robust_se(U:np.ndarray, V:np.ndarray, lamb:float=100, t_cutoff:float=6.5, abs_t:bool=False) -> sparse.csc_matrix:
    from . import _core as core
    if U.dtype == np.float32 and V.dtype == np.float32:
        return core.py_robust_se_ff(U, V, lamb, t_cutoff, abs_t)
    elif U.dtype == np.float64 and V.dtype == np.float64:
        return core.py_robust_se_dd(U, V, lamb, t_cutoff, abs_t)
    elif U.dtype == np.float32 and V.dtype == np.float64:
        return core.py_robust_se_fd(U, V, lamb, t_cutoff, abs_t)
    elif U.dtype == np.float64 and V.dtype == np.float32:
        return core.py_robust_se_df(U, V, lamb, t_cutoff, abs_t)
    else:
        raise NotImplementedError("Not implemented for numpy dtype")

def robust_prepare(U:np.ndarray, V:np.ndarray, B=None, n_components:int=None, min_norm:float=1e-5, return_U=True, power:float=0):
    if U.shape[1] != V.shape[0]:
        raise ValueError("Shapes must match")
    if n_components != None:
        if n_components <= 0:
            raise ValueError("N components must be > 0")
        U = U[:, :n_components]
        V = V[:n_components, :]
    bad_flag = np.linalg.norm(V, axis=0, ord=2) < min_norm
    if B is None:
        return U, V[:, ~bad_flag], ~bad_flag
    if U.shape[0] != B.shape[0]:
        raise ValueError("U and B are not the correct shape")
    V1 = np.nan * np.ones_like(V)
    V_svd = np.linalg.svd(V[:, ~bad_flag], full_matrices=False)
    U = U @ V_svd.U
    qr = np.linalg.qr(U - B @ (np.linalg.pinv(B) @ U))
    RS_svd = np.linalg.svd(qr.R @ np.diag(V_svd.S))
    V1[:, ~bad_flag] = np.diag(RS_svd.S**power) @ RS_svd.Vh @ V_svd.Vh
    if return_U:
        return qr.Q @ RS_svd.U @ np.diag(RS_svd.S**(1-power)), V1, ~bad_flag
    else:
        return None, V1, ~bad_flag

def robust_se_default(U, V, B=None, t_cutoff:float=None, abs_t:bool=False, lamb:float=100., nominal_p_cutoff:float=0.05, n_components:int=None, min_norm:float=1e-5) -> sparse.csc_matrix:
    """
    U: U from SVD
    V: V\Sigma from SVD
    B: Batch matrix from pd.get_dummies, or just intercept.
    """
    if U.shape[1] != V.shape[0]:
        raise ValueError("U and V must have compatible dimensions : %d != %d" % (U.shape[1], V.shape[0]))
    U, V, var_flag = robust_prepare(U, V, B=B, n_components=n_components, min_norm=min_norm, return_U=True)
    if t_cutoff is None:
        import scipy.stats
        t_cutoff = scipy.stats.t.isf(nominal_p_cutoff * V.shape[1]**-2,
                                     max(U.shape[0] - V.shape[1] - 2, 1))
    S = _robust_se(U, V, t_cutoff=t_cutoff, abs_t=abs_t, lamb=lamb)
    I = np.ravel(np.where(var_flag))
    trans = sparse.csr_matrix((np.ones_like(I), (I, np.arange(len(I)))),
                              shape=(len(var_flag), S.shape[0]))
    return trans.dot(S).dot(trans.T)


