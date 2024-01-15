import scdemon_ext

def robust_se(U, V, B=None, t_cutoff:float=None, abs_t:bool=False, nominal_p_cutoff:float=0.05):
    """
    U: U from SVD
    V: V\Sigma from SVD
    B: Batch matrix from pd.get_dummies, or just intercept.
    """
    import numpy as np
    if U.shape[1] != V.shape[0]:
        raise ValueError("U and V must have compatible dimensions : %d != %d" % (U.shape[1], V.shape[0]))
    if B is None:
        B = np.ones((U.shape[0], 1), dtype="f8")
    elif U.shape[0] != B.shape[0]:
        raise ValueError("U and B must have compatible dimensions : %d != %d" % (U.shape[0], B.shape[0]))
    UpB = scdemon_ext.ols_beta(U.astype("f8"), B.astype("f8"))
    UpU = scdemon_ext.ols_beta(U.astype("f8"), U.astype("f8"))
    if t_cutoff is None:
        import scipy.stats
        t_cutoff = scipy.stats.t.isf(nominal_p_cutoff * V.shape[1]**-2,
                                     B.shape[0] - B.shape[1])
    M = scdemon_ext.robust_se(V.astype("f8"), UpU, UpB, t_cutoff, abs_t)
    return M

    
