import numpy as np
import scipy.sparse

def _robust_se(V:np.ndarray, t_cutoff:float, abs_t:bool) -> scipy.sparse.csc_matrix:
    from ._core import py_robust_se
    from .progress_manager import ProgressManager, _interrupt_checker
    pm = ProgressManager()
    return py_robust_se(V.astype('f8'), pm.hook, _interrupt_checker, t_cutoff, abs_t)
    
def _robust_se_X(i:int, V:np.ndarray) -> np.ndarray:
    from ._core import py_robust_se_X
    return py_robust_se_X(i, V.astype("f8"))

def robust_se(U, V, B=None, t_cutoff:float=None, abs_t:bool=False, nominal_p_cutoff:float=0.05) -> scipy.sparse.csc_matrix:
    """
    U: U from SVD
    V: V\Sigma from SVD
    B: Batch matrix from pd.get_dummies, or just intercept.
    """
    
    from ._core import ols_beta
    if U.shape[1] != V.shape[0]:
        raise ValueError("U and V must have compatible dimensions : %d != %d" % (U.shape[1], V.shape[0]))
    if B is None:
        B = np.ones((U.shape[0], 1), dtype="f8")
    elif U.shape[0] != B.shape[0]:
        raise ValueError("U and B must have compatible dimensions : %d != %d" % (U.shape[0], B.shape[0]))
    UpB = ols_beta(U.astype("f8"), B.astype("f8"))
    UpU = ols_beta(U.astype("f8"), U.astype("f8"))
    if t_cutoff is None:
        import scipy.stats
        t_cutoff = scipy.stats.t.isf(nominal_p_cutoff * V.shape[1]**-2,
                                     B.shape[0] - B.shape[1])
    return _robust_se(V, t_cutoff=t_cutoff, abs_t=abs_t)


