# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import division, print_function

from math import pi

import numpy as np
import tensorflow as tf

from .pymbd import freq_grid

__all__ = ()

pi = tf.constant(pi, tf.float64)
ang = 1 / 0.529177249


class MBDEvaluator(object):
    def __init__(self, gradients=False, **kwargs):
        self._inputs = coords, alpha_0, C6, R_vdw, beta = [
            tf.placeholder(tf.float64, shape=shape, name=name)
            for shape, name in [
                ((None, 3), 'coords'),
                ((None,), 'alpha_0'),
                ((None,), 'C6'),
                ((None,), 'R_vdw'),
                ((), 'beta'),
            ]
        ]
        self._output = mbd_energy(*self._inputs, **kwargs)
        if gradients:
            self._init_gradients()
        else:
            self._gradients = None

    def _init_gradients(self):
        self._gradients = tf.gradients(self._output, [self._inputs[0]])[0]

    def __call__(self, coords, alpha_0, C6, R_vdw, beta=0.83, gradients=None):
        inputs = dict(zip(self._inputs, [coords, alpha_0, C6, R_vdw, beta]))
        outputs = self._output
        if gradients or gradients is None and self._gradients is not None:
            if self._gradients is None:
                self._init_gradients()
            outputs = outputs, self._gradients
        return tf.get_default_session().run(outputs, inputs)


def mbd_energy(coords, alpha_0, C6, R_vdw, beta, nfreq=15):
    freq, freq_w = freq_grid(nfreq)
    omega = 4 / 3 * C6 / alpha_0 ** 2
    alpha_dyn = [alpha_0 / (1 + (u / omega) ** 2) for u in freq]
    alpha_dyn_rsscs = []
    for a in alpha_dyn:
        sigma = (tf.sqrt(2 / pi) * a / 3) ** (1 / 3)
        dipmat = dipole_matrix(
            coords, 'fermi,dip,gg', sigma=sigma, R_vdw=R_vdw, beta=beta
        )
        a_nlc = tf.linalg.inv(tf.diag(_repeat(1 / a, 3)) + dipmat)
        a_contr = sum(tf.reduce_sum(a_nlc[i::3, i::3], 1) for i in range(3)) / 3
        alpha_dyn_rsscs.append(a_contr)
    alpha_dyn_rsscs = tf.stack(alpha_dyn_rsscs)
    C6_rsscs = 3 / pi * tf.reduce_sum(freq_w[:, None] * alpha_dyn_rsscs ** 2, 0)
    R_vdw_rsscs = R_vdw * (alpha_dyn_rsscs[0, :] / alpha_0) ** (1 / 3)
    omega_rsscs = 4 / 3 * C6_rsscs / alpha_dyn_rsscs[0, :] ** 2
    dipmat = dipole_matrix(coords, 'fermi,dip', R_vdw=R_vdw_rsscs, beta=beta)
    pre = _repeat(omega_rsscs * tf.sqrt(alpha_dyn_rsscs[0, :]), 3)
    eigs = tf.linalg.eigvalsh(
        tf.diag(_repeat(omega_rsscs ** 2, 3)) + pre[:, None] * pre[None, :] * dipmat
    )
    ene = tf.reduce_sum(tf.sqrt(eigs)) / 2 - 3 * tf.reduce_sum(omega_rsscs) / 2
    return ene


def dipole_matrix(coords, damping, beta=0.0, R_vdw=None, sigma=None, a=6.0):
    Rs = coords[:, None, :] - coords[None, :, :]
    if R_vdw is not None:
        S_vdw = beta * (R_vdw[:, None] + R_vdw[None, :])
        # 1e10 rather than inf is necessary to avoid nan gradients
        dists = tf.sqrt(_set_diag(tf.reduce_sum(Rs ** 2, -1), 1e10))
    if sigma is not None:
        sigmaij = tf.sqrt(sigma[:, None] ** 2 + sigma[None, :] ** 2)
    if damping == 'fermi,dip':
        dipmat = damping_fermi(dists, S_vdw, a)[:, :, None, None] * T_bare(Rs)
    elif damping == 'fermi,dip,gg':
        dipmat = (1 - damping_fermi(dists, S_vdw, a)[:, :, None, None]) * T_erf_coulomb(
            Rs, sigmaij
        )
    else:
        raise ValueError(f'Unsupported damping: {damping}')
    n_atoms = tf.shape(coords)[0]
    return tf.reshape(tf.transpose(dipmat, (0, 2, 1, 3)), (3 * n_atoms, 3 * n_atoms))


def damping_fermi(R, S_vdw, d):
    return 1 / (1 + tf.exp(-d * (R / S_vdw - 1)))


def T_bare(R):
    R_2 = tf.reduce_sum(R ** 2, -1)
    # 1e10 rather than inf is necessary to avoid nan gradients
    R_5 = tf.sqrt(_set_diag(R_2, 1e10)) ** 5
    return (
        -3 * R[:, :, :, None] * R[:, :, None, :]
        + R_2[:, :, None, None] * np.eye(3)[None, None, :, :]
    ) / R_5[:, :, None, None]


def T_erf_coulomb(R, sigma):
    bare = T_bare(R)
    # 1e-10 rather than 0 is necessary to avoid nan gradients from sqrt(0)
    R_1 = tf.sqrt(_set_diag(tf.reduce_sum(R ** 2, -1), 1e-10))
    # 1e10 rather than inf is necessary to avoid nan gradients
    R_5 = _set_diag(R_1 ** 5, 1e10)
    RR_R5 = R[:, :, :, None] * R[:, :, None, :] / R_5[:, :, None, None]
    zeta = R_1 / sigma
    theta = 2 * zeta / tf.sqrt(pi) * tf.exp(-(zeta ** 2))
    erf_theta = tf.erf(zeta) - theta
    return (
        erf_theta[:, :, None, None] * bare
        + (2 * (zeta ** 2) * theta)[:, :, None, None] * RR_R5
    )


def _set_diag(A, val):
    return tf.matrix_set_diag(A, tf.fill(tf.shape(A)[0:1], tf.cast(val, tf.float64)))


def _repeat(a, n):
    return tf.reshape(tf.tile(a[:, None], (1, n)), (-1,))
