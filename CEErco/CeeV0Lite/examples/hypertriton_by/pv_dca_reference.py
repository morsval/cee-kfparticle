"""Formal default reference preparation, independent of MC truth."""
from pathlib import Path
import ctypes
import numpy as np

def prepare_reference(x,c,pv,by=.5):
    lib=ctypes.CDLL(str(Path(__file__).resolve().parent/'libtoy_core.so'))
    P=ctypes.POINTER(ctypes.c_double)
    lib.prepare_reference_batch.argtypes=[ctypes.c_int,P,P,P,ctypes.c_double,P,P,P]
    args=[np.ascontiguousarray(v,dtype=np.float64) for v in (x,c,pv)]
    xo=np.empty_like(args[0]);co=np.empty_like(args[1]);diag=np.empty((len(x),2,4))
    lib.prepare_reference_batch(len(x),*[v.ctypes.data_as(P) for v in args],by,
        xo.ctypes.data_as(P),co.ctypes.data_as(P),diag.ctypes.data_as(P))
    return xo,co,diag
