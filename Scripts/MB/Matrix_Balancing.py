import sys

import numpy

def BNEWT(A, tol=1e-6, x0=None, delta=0.1, Delta=3, fl=False):
    # Initialize values
    res = []
    n = A.shape[0]
    e = numpy.ones((n, 1), dtype=numpy.float64)
    if x0 is None:
        x0=numpy.copy(e)
    g = 0.9
    eta = etamax = 0.1
    stop_tol = tol * 0.5
    x = numpy.copy(x0)
    rt = tol ** 2.0
    v = x * numpy.dot(A, x)
    rk = 1.0 - v
    rho_km1 = numpy.dot(rk.T, rk)[0, 0]
    rho_km2 = rho_km1
    rold = rout = rho_km1
    i = MVP = 0
    if fl:
        print >> sys.stderr, ("it in. it res stop cur\n"),
    # Outer iteration
    while rout > rt:
        i += 1
        k = 0
        y = numpy.copy(e)
        innertol = max(eta ** 2.0 * rout, rt)
        # Inner iteration by CG
        while rho_km1 > innertol:
            k += 1
            if k == 1:
                Z = rk / v
                p = numpy.copy(Z)
                rho_km1 = numpy.dot(rk.T, Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p
            # Update search direction efficiently
            w = x * numpy.dot(A, x * p) + v * p
            alpha = rho_km1 / numpy.dot(p.T, w)[0, 0]
            ap = alpha * p
            # Test distance to boundary of cone
            ynew = y + ap
            if numpy.amin(ynew) <= delta:
                if delta == 0:
                    break
                ind = numpy.where(ap < 0.0)[0]
                gamma = numpy.amin((delta - y[ind]) / ap[ind])
                y += gamma * ap
                break
            if numpy.amax(ynew) >= Delta:
                ind = numpy.where(ynew > Delta)[0]
                gamma = numpy.amin((Delta - y[ind]) / ap[ind])
                y += gamma * ap
                break
            y = numpy.copy(ynew)
            rk -= alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            rho_km1 = numpy.dot(rk.T, Z)[0, 0]
        x *= y
        v = x * numpy.dot(A, x)
        rk = 1.0 - v
        rho_km1 = numpy.dot(rk.T, rk)[0, 0]
        rout = rho_km1
        MVP += k + 1
        # Update inner iteration stopping criterion
        rat = rout / rold
        rold = rout
        res_norm = rout ** 0.5
        eta_o = eta
        eta = g * rat
        if g * eta_o ** 2.0 > 0.1:
            eta = max(eta, g * eta_o ** 2.0)
        eta = max(min(eta, etamax), stop_tol / res_norm)
        if fl:
            print >> sys.stderr, ("%03i %06i %03.3f %e %e\n") % (i, k, res_norm, rt, rout),
            res.append(res_norm)
    if fl:
        print >> sys.stderr, ("Matrix-vector products = %06i\n") % (MVP),
    return [x, res]


def BNEWT_sparse(A, tol=1e-6, x0=None, delta=0.1, Delta=3, fl=False):
    # A should be an n by 3 array of index1, index2, count
    # Initialize values
    n = numpy.amax(A[:, :2]) + 1
    res = []
    e = numpy.ones((n, 1), dtype=numpy.float64)
    if x0 is None:
        x0=numpy.copy(e)
    g = 0.9
    eta = etamax = 0.1
    stop_tol = tol * 0.5
    x = numpy.copy(x0)
    rt = tol ** 2.0
    corr = (A[:, 2:3] * x[A[:, 0]] * x[A[:, 1]]).reshape(-1)
    v = numpy.bincount(A[:, 0], weights=corr, minlength=n).reshape(-1, 1)
    v += numpy.bincount(A[:, 1], weights=corr, minlength=n).reshape(-1, 1)
    rk = 1.0 - v
    rho_km1 = numpy.dot(rk.T, rk)[0, 0]
    rho_km2 = rho_km1
    rold = rout = rho_km1
    i = MVP = 0
    if fl:
        print >> sys.stderr, ("it in. it res stop cur\n"),
    # Outer iteration
    while rout > rt:
        i += 1
        k = 0
        y = numpy.copy(e)
        innertol = max(eta ** 2.0 * rout, rt)
        # Inner iteration by CG
        while rho_km1 > innertol:
            k += 1
            if k == 1:
                Z = rk / v
                p = numpy.copy(Z)
                rho_km1 = numpy.dot(rk.T, Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p
            # Update search direction efficiently
            corr = A[:, 2:3] * x[A[:, 0]] * x[A[:, 1]]
            w = numpy.bincount(A[:, 0], weights=(corr * p[A[:, 1]]).reshape(-1), minlength=n).reshape(-1, 1)
            w += numpy.bincount(A[:, 1], weights=(corr * p[A[:, 0]]).reshape(-1), minlength=n).reshape(-1, 1)
            w += v * p
            alpha = rho_km1 / numpy.dot(p.T, w)[0, 0]
            ap = alpha * p
            # Test distance to boundary of cone
            ynew = y + ap
            if numpy.amin(ynew) <= delta:
                if delta == 0:
                    break
                ind = numpy.where(ap < 0.0)[0]
                gamma = numpy.amin((delta - y[ind]) / ap[ind])
                y += gamma * ap
                break
            if numpy.amax(ynew) >= Delta:
                ind = numpy.where(ynew > Delta)[0]
                gamma = numpy.amin((Delta - y[ind]) / ap[ind])
                y += gamma * ap
                break
            y = numpy.copy(ynew)
            rk -= alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            rho_km1 = numpy.dot(rk.T, Z)[0, 0]
        x *= y
        corr = (A[:, 2:3] * x[A[:, 0]] * x[A[:, 1]]).reshape(-1)
        v = numpy.bincount(A[:, 0], weights=corr, minlength=n).reshape(-1, 1)
        v += numpy.bincount(A[:, 1], weights=corr, minlength=n).reshape(-1, 1)
        rk = 1.0 - v
        rho_km1 = numpy.dot(rk.T, rk)[0, 0]
        rout = rho_km1
        MVP += k + 1
        # Update inner iteration stopping criterion
        rat = rout / rold
        rold = rout
        res_norm = rout ** 0.5
        eta_o = eta
        eta = g * rat
        if g * eta_o ** 2.0 > 0.1:
            eta = max(eta, g * eta_o ** 2.0)
        eta = max(min(eta, etamax), stop_tol / res_norm)
        if fl:
            print >> sys.stderr, ("%03i %06i %03.3f %e %e\n") % (i, k, res_norm, rt, rout),
            res.append(res_norm)
    if fl:
        print >> sys.stderr, ("Matrix-vector products = %06i\n") % (MVP),
    return [x, res]


def BNEWT_sparse_binary(A, tol=1e-6, x0=None, delta=0.1, Delta=3, fl=False):
    # A should be an n by 3 array of index1, index2, count
    # Initialize values
    n = numpy.amax(A[:, :2]) + 1
    res = []
    e = numpy.ones((n, 1), dtype=numpy.float64)
    if x0 is None:
        x0=numpy.copy(e)
    g = 0.9
    eta = etamax = 0.1
    stop_tol = tol * 0.5
    x = numpy.copy(x0)
    rt = tol ** 2.0
    corr = (x[A[:, 0]] * x[A[:, 1]]).reshape(-1)
    v = numpy.bincount(A[:, 0], weights=corr, minlength=n).reshape(-1, 1)
    v += numpy.bincount(A[:, 1], weights=corr, minlength=n).reshape(-1, 1)
    rk = 1.0 - v
    rho_km1 = numpy.dot(rk.T, rk)[0, 0]
    rho_km2 = rho_km1
    rold = rout = rho_km1
    i = MVP = 0
    if fl:
        print >> sys.stderr, ("it in. it res stop cur\n"),
    # Outer iteration
    while rout > rt:
        i += 1
        k = 0
        y = numpy.copy(e)
        innertol = max(eta ** 2.0 * rout, rt)
        # Inner iteration by CG
        while rho_km1 > innertol:
            k += 1
            if k == 1:
                Z = rk / v
                p = numpy.copy(Z)
                rho_km1 = numpy.dot(rk.T, Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p
            # Update search direction efficiently
            corr = x[A[:, 0]] * x[A[:, 1]]
            w = numpy.bincount(A[:, 0], weights=(corr * p[A[:, 1]]).reshape(-1), minlength=n).reshape(-1, 1)
            w += numpy.bincount(A[:, 1], weights=(corr * p[A[:, 0]]).reshape(-1), minlength=n).reshape(-1, 1)
            w += v * p
            alpha = rho_km1 / numpy.dot(p.T, w)[0, 0]
            ap = alpha * p
            # Test distance to boundary of cone
            ynew = y + ap
            if numpy.amin(ynew) <= delta:
                if delta == 0:
                    break
                ind = numpy.where(ap < 0.0)[0]
                gamma = numpy.amin((delta - y[ind]) / ap[ind])
                y += gamma * ap
                break
            if numpy.amax(ynew) >= Delta:
                ind = numpy.where(ynew > Delta)[0]
                gamma = numpy.amin((Delta - y[ind]) / ap[ind])
                y += gamma * ap
                break
            y = numpy.copy(ynew)
            rk -= alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            rho_km1 = numpy.dot(rk.T, Z)[0, 0]
        x *= y
        corr = (x[A[:, 0]] * x[A[:, 1]]).reshape(-1)
        v = numpy.bincount(A[:, 0], weights=corr, minlength=n).reshape(-1, 1)
        v += numpy.bincount(A[:, 1], weights=corr, minlength=n).reshape(-1, 1)
        rk = 1.0 - v
        rho_km1 = numpy.dot(rk.T, rk)[0, 0]
        rout = rho_km1
        MVP += k + 1
        # Update inner iteration stopping criterion
        rat = rout / rold
        rold = rout
        res_norm = rout ** 0.5
        eta_o = eta
        eta = g * rat
        if g * eta_o ** 2.0 > 0.1:
            eta = max(eta, g * eta_o ** 2.0)
        eta = max(min(eta, etamax), stop_tol / res_norm)
        if fl:
            print >> sys.stderr, ("%03i %06i %03.3f %e %e\n") % (i, k, res_norm, rt, rout),
            res.append(res_norm)
    if fl:
        print >> sys.stderr, ("Matrix-vector products = %06i\n") % (MVP),
    return [x, res]

