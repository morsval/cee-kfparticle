#!/usr/bin/env python3
"""Compare every reconstructed toy candidate to the previous Python fit.
The reference uses a 3x3 inverse with weighted elimination of S, whereas the
formal C++ core uses a 2D transverse projection and a square-root transport.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import json
from pathlib import Path
import numpy as np
from reference_charged_topology import fit as reference_fit
from simulate import PV_SIGMA, BY

root = Path(__file__).resolve().parent
data = np.load(root/'output/candidates.npz')
indices = np.flatnonzero(data['base'])
counts = {'compared':0, 'both_valid':0, 'cpp_only':0, 'python_only':0, 'both_invalid':0}
errors = {'state_in_sigma':0., 'covariance_scaled':0., 'chi2':0.}
for begin in range(0, len(indices), 5000):
    ii = indices[begin:begin+5000]
    q = data['raw'][ii, 3:10]
    c = data['raw'][ii, 10:59].reshape(-1, 7, 7)
    pvc = np.broadcast_to(np.eye(3)*PV_SIGMA**2, (len(ii), 3, 3))
    ref, refc, refchi, ok, _ = reference_fit(q, c, data['pv'][ii], pvc, BY)
    cppok = data['topo_valid'][ii]
    same = ok & cppok
    counts['compared'] += len(ii)
    counts['both_valid'] += int(same.sum())
    counts['cpp_only'] += int((cppok & ~ok).sum())
    counts['python_only'] += int((ok & ~cppok).sum())
    counts['both_invalid'] += int((~ok & ~cppok).sum())
    if same.any():
        scale = np.sqrt(np.diagonal(refc[same], axis1=1, axis2=2))
        errors['state_in_sigma'] = max(errors['state_in_sigma'], float(np.max(abs(data['q8'][ii[same]]-ref[same])/scale)))
        errors['covariance_scaled'] = max(errors['covariance_scaled'], float(np.max(abs(data['c8'][ii[same]]-refc[same])/(scale[:, :, None]*scale[:, None, :]))))
        errors['chi2'] = max(errors['chi2'], float(np.max(abs(data['chi'][ii[same]]-refchi[same]))))

assert counts['python_only'] == 0, counts
assert errors['state_in_sigma'] < 2e-5, errors
assert errors['covariance_scaled'] < 2e-5, errors
assert errors['chi2'] < 1e-5, errors
result = dict(counts=counts, maximum_differences=errors, passed=True)
(root/'output/formal_vs_python.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
