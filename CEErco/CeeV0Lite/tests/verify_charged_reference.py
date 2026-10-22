#!/usr/bin/env python3
"""Independent charged validation; NumPy/SciPy only, no production Python code.

Eight production-end parameters are fitted to ten correlated measurements.
Numerical model derivatives and a Lorentz-force ODE provide separate checks
of the C++ constraint algebra, covariance and magnetic-field signs.
"""
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

root = Path(__file__).resolve().parents[1]
raw = subprocess.check_output([str(root/'obj/test_charged_production'), '--reference'], text=True)
cases = [json.loads(line) for line in raw.splitlines()]
max_state = max_cov = max_chi2 = max_ode = 0.0
for data in cases:
    y = np.array(data['input'])
    c = np.array(data['covariance']).reshape(10, 10)
    chol = np.linalg.cholesky(c)
    charge, by = data['charge'], data['by']
    curvature = 0.00299792458*charge*by

    def model(q):
        s = q[7]
        angle = curvature*s
        si = s*np.sinc(angle/np.pi)
        om = s*(angle/2)*np.sinc(angle/(2*np.pi))**2
        px, py, pz = q[3:6]
        decay_position = q[:3]+np.array([si*px-om*pz, s*py, si*pz+om*px])
        decay_momentum = np.array([np.cos(angle)*px-np.sin(angle)*pz,
                                   py, np.sin(angle)*px+np.cos(angle)*pz])
        return np.r_[decay_position, decay_momentum, q[6], q[:3]]

    s0 = np.dot(y[:3]-y[7:10], y[3:6])/np.dot(y[3:6], y[3:6])
    a0 = curvature*s0
    p0 = np.array([np.cos(a0)*y[3]+np.sin(a0)*y[5], y[4],
                   -np.sin(a0)*y[3]+np.cos(a0)*y[5]])
    start = np.r_[y[7:10], p0, y[6], s0]
    fit = least_squares(lambda q: solve_triangular(chol, model(q)-y, lower=True),
                        start, jac='3-point', x_scale='jac', xtol=1e-13,
                        ftol=1e-13, gtol=1e-12, max_nfev=300)
    assert fit.success, fit.message
    ref_cov = np.linalg.inv(fit.jac.T@fit.jac)
    cpp = np.array(data['state'])
    cpp_cov = np.array(data['state_covariance']).reshape(8, 8)
    scale = np.sqrt(np.diag(ref_cov))
    max_state = max(max_state, float(np.max(abs(fit.x-cpp)/scale)))
    max_cov = max(max_cov, float(np.max(abs(ref_cov-cpp_cov)/np.outer(scale, scale))))
    max_chi2 = max(max_chi2, abs(2*fit.cost-data['chi2']))

    def rhs(length, state):
        momentum = state[3:6]
        magnitude = np.linalg.norm(momentum)
        return np.r_[momentum/magnitude,
                     0.00299792458*charge*np.cross(momentum, [0., by, 0.])/magnitude]

    length = cpp[7]*np.linalg.norm(cpp[3:6])
    integrated = solve_ivp(rhs, (0., length), cpp[:6], rtol=2e-12, atol=1e-13).y[:, -1]
    max_ode = max(max_ode, float(np.max(abs(integrated-np.array(data['decay_state'])[:6]))))

assert len(cases) == 16
assert max_state < 2e-5, max_state
assert max_cov < 2e-6, max_cov
assert max_chi2 < 1e-8, max_chi2
assert max_ode < 2e-8, max_ode
print(json.dumps({'cases':len(cases), 'max_state_difference_in_sigma':max_state,
                  'max_scaled_covariance_difference':max_cov,
                  'max_chi2_difference':max_chi2,
                  'max_Lorentz_ODE_difference':max_ode, 'passed':True}, indent=2))
