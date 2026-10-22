"""Thin Python adapter: all production-vertex fits run the formal C++ core.

This example specifies the hypertriton charge +1 explicitly. The independent
Python implementation is retained only as a numerical reference and for
plotting derived quantities, never as a reconstruction fallback.
"""
import ctypes
from pathlib import Path
import numpy as np
from reference_charged_topology import KAPPA, helix_matrices, observables, endpoint

def fit(q, c, pv, pvc, by=.5, maxiter=30):
    lib = ctypes.CDLL(str(Path(__file__).resolve().parent/'libtoy_core.so'))
    pointer = ctypes.POINTER(ctypes.c_double)
    lib.production_fit_batch.argtypes = [ctypes.c_int, pointer, pointer, pointer,
                                          pointer, ctypes.c_int, ctypes.c_double,
                                          ctypes.c_int, pointer]
    lib.production_fit_batch.restype = None
    q, c, pv, pvc = [np.ascontiguousarray(x, dtype=np.float64) for x in (q, c, pv, pvc)]
    n = len(q)
    assert q.shape == (n, 7) and c.shape == (n, 7, 7)
    assert pv.shape == (n, 3) and pvc.shape == (n, 3, 3)
    out = np.empty((n, 75))
    lib.production_fit_batch(n, q.ctypes.data_as(pointer), c.ctypes.data_as(pointer),
                             pv.ctypes.data_as(pointer), pvc.ctypes.data_as(pointer),
                             1, by, maxiter, out.ctypes.data_as(pointer))
    return out[:, :8], out[:, 8:72].reshape(n, 8, 8), out[:, 72], out[:, 73] == 1, out[:, 74].astype(int)
