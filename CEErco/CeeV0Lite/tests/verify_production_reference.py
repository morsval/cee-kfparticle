#!/usr/bin/env python3
"""Optional independent cross-check, requiring NumPy/SciPy.

Fit eight physical parameters directly to ten correlated measurements,
using SciPy least_squares rather than the production code's 10D constraint.
"""
import json
import pathlib
import subprocess
import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

root = pathlib.Path(__file__).resolve().parents[1]
data = json.loads(subprocess.check_output([str(root/'obj/test_production_state'), '--reference'], text=True))
y = np.array(data['input'])
c = np.array(data['covariance']).reshape(10, 10)
chol = np.linalg.cholesky(c)

def model(q):
    return np.r_[q[:3]+q[7]*q[3:6], q[3:7], q[:3]]

def jac(q):
    h = np.zeros((10, 8))
    h[:3,:3] = np.eye(3)
    h[:3,3:6] = q[7]*np.eye(3)
    h[:3,7] = q[3:6]
    h[3:7,3:7] = np.eye(4)
    h[7:10,:3] = np.eye(3)
    return h

start = np.r_[y[7:10], y[3:7], np.dot(y[:3]-y[7:10],y[3:6])/np.dot(y[3:6],y[3:6])]
fit = least_squares(lambda q: solve_triangular(chol, model(q)-y, lower=True), start,
                    jac=lambda q: solve_triangular(chol, jac(q), lower=True),
                    xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=300)
assert fit.success, fit.message
np.testing.assert_allclose(data['state'], fit.x, rtol=2e-7, atol=2e-8)
reference_covariance = np.linalg.inv(fit.jac.T@fit.jac)
np.testing.assert_allclose(np.array(data['state_covariance']).reshape(8,8),
                           reference_covariance, rtol=2e-5, atol=2e-9)
np.testing.assert_allclose(data['chi2'], 2*fit.cost, rtol=2e-7, atol=1e-9)
print('Independent 8-parameter GLS fit agrees: state, full covariance and chi2')
