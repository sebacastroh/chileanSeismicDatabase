# -*- coding: utf-8 -*-
"""
Created on Sat May 15 16:51:41 2021

@author: sebac
"""
import cython
from libc.math cimport fabs, exp, log, sqrt, pow
from libc.stdio cimport printf, fflush, stdout

cdef extern from "<math.h>":
    const double M_PI
    const double INFINITY

@cython.boundscheck(False)
@cython.wraparound(False)
cdef long choose_interval(double[::1] cdf, long n, double u):
    cdef long i
    
    for i in range(n):
        if u < cdf[i]:
            return i
        
    return -1

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void mat_zero(double[:,::1] A, long n):
    """
    Converts matrix to a zero matrix
    
    Parameters
    ----------
    A : double**
        Matrix to be converted
    n : long
        Size of the matrix
        
    Returns
    -------
    None - Matrix A is modified
    """
    cdef long i
    cdef long j
    
    for i in range(n):
        for j in range(n):
            A[i][j] = 0.

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double mat_pivot(double[:,::1] A, double[:,::1] P, long n):
    """
    Computes the pivoting matrix of an A matrix
    
    Parameters
    ----------
    A : double**
        Matrix to be pivoted
    P : double**
        Pivot matrix
    n : long
        Size of the matrix
        
    Returns
    -------
    p : double
        Number of permutations    
    
    Matrix P is modified
    """
    cdef long i
    cdef long j
    cdef long k
    cdef long max_j
    cdef double temp
    
    cdef double p
    
    for i in range(n):
        for j in range(n):
            P[i][j] = (i == j)
            
    p = 0.
    
    for i in range(n):
        max_j = i
        for j in range(i, n):
            if fabs(A[j][i]) > fabs(A[max_j][i]):
                max_j = j
        
        if max_j != i:
            p += 1.
            for k in range(n):
                temp = P[i][k]
                P[i][k] = P[max_j][k]
                P[max_j][k] = temp
    
    return p

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void mat_mul(double[:,::1] A, double[:,::1] B, long n, double[:,::1] C):
    """
    Performs A*B multiplication. It assumes that A and B are square matrices of
    the same size
    
    Parameters
    ----------
    A : double**
        Left matrix of the multiplication
    B : double**
        Right matrix of the multiplication
    n : long
        Size of the matrices
    C : double**
        Resulting matrix
        
    Returns
    -------
    None - Matrix C is modified
    """
    cdef long i
    cdef long j
    cdef long k
    
    for i in range(n):
        for j in range(n):
            C[i][j] = 0.
            for k in range(n):
                C[i][j] += A[i][k]*B[k][j]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double mat_LU(double[:,::1] A, double[:,::1] L, double[:,::1] U, double[:,::1] P, long n, double[:,::1] Aprime):
    """
    Performs P*A = L*U factorization given an A matrix of size n
    
    Parameters
    ----------
    A : double**
        Original matrix
    L : double**
        Lower triangular factorization matrix
    U : double**
        Upper triangular factorization matrix
    P : double**
        Pivot factorization matrix
    n : long
        Size of the matrices
    Aprime : double**
        Auxiliar matrix
        
    Returns
    -------
    det : double
        Determinant of A
    
    Matrix C is modified
    """
    cdef long i
    cdef long j
    cdef long k
    cdef double s
    cdef double p
    cdef double det
    
    mat_zero(L, n)
    mat_zero(U, n)
    mat_zero(P, n)
    p = mat_pivot(A, P, n)
    
    mat_mul(P, A, n, Aprime)
    
    for i in range(n):
        L[i][i] = 1.
        
    det = 1.
    for i in range(n):
        for j in range(n):
            if j <= i:
                s = 0.
                for k in range(j):
                    s += L[j][k]*U[k][i]
                U[j][i] = Aprime[j][i] - s
            
            if j >= i:
                s = 0.
                for k in range(i):
                    s += L[j][k]*U[k][i]
                L[j][i] = (Aprime[j][i] - s)/U[i][i]
        det *= U[i][i]
    
    det *= pow(-1, p)
    
    return det

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void mat_inverse_triangular_lower(double[:,::1] L, double[:,::1] Li, long n):
    """
    Computes the inverse of a lower triangular matrix
    
    Parameters
    ----------
    L : double**
        Lower triangular matrix
    Li : double**
        Inverse lower triangular matrix
    n : long
        Size of the matrices
        
    Returns
    -------
    None - Li is modified
    """
    cdef long i
    cdef long j
    cdef long k
    
    mat_zero(Li, n)
    
    for k in range(n):
        Li[k][k] = 1./L[k][k]
        for i in range(k+1, n):
            Li[i][k] = 0
            for j in range(k, i):
                Li[i][k] -= L[i][j]*Li[j][k]
            Li[i][k] /= L[i][i]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void mat_inverse_triangular_upper(double[:,::1] U, double[:,::1] Ui, long n):
    """
    Computes the inverse of an upper triangular matrix
    
    Parameters
    ----------
    U : double**
        Upper triangular matrix
    Ui : double**
        Inverse upper triangular matrix
    n : long
        Size of the matrices
        
    Returns
    -------
    None - Ui is modified
    """
    cdef long i
    cdef long j
    cdef long k
    
    mat_zero(Ui, n)
    
    for k in range(n-1, -1, -1):
        Ui[k][k] = 1./U[k][k]
        for i in range(k-1, -1, -1):
            Ui[i][k] = 0.
            for j in range(k, i, -1):
                Ui[i][k] -= U[i][j]*Ui[j][k]
            Ui[i][k] /= U[i][i]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef long Cholesky_Decomposition(double[:,::1] S, double[:,::1] B, long n):
    """
    Computes the Cholesky decomposition of a S matrix, i.e. S = B*B^T
    
    Parameters
    ----------
    S : double**
        Original matrix
    B : double**
        Cholesky decomposition
    n : long
        Size of the matrices
        
    Returns
    -------
    valid : long
        If valid is equal to 1, the Cholesky was performed successfully, in other
        case the original matrix is not a positive-definite matrix
    
    B is modified
    """
    cdef long i
    cdef long j
    cdef long k
    
    cdef double s
    cdef double det
    
    mat_zero(B, n)
    
    for i in range(n):
        for j in range(i+1):
            s = 0.
            if i == j:
                for k in range(j):
                    s += B[j][k]*B[j][k]
                B[j][j] = sqrt(S[j][j] - s)
            else:
                for k in range(j):
                    s += B[i][k]*B[j][k]
                if B[j][j] > 0.:
                    B[i][j] = (S[i][j] - s)/B[j][j]
                else:
                    return -1
                
    return 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double update_partition(double[:,::1] data, double[::1] fg, double[::1] fl, long xi, long xj, long kmax,
                              long[::1] nlocal_parameters, double[:,::1] mu, double[:,::1] sigma, double[:,:,::1] S,
                              double[::1] detS, double[::1] epsilon, double[:,:,::1] X, double[:,::1] N, double[:,::1] Q,
                              double[:,::1] L, double[:,::1] U, double[:,::1] P, double[:,::1] aux1, double[:,::1] aux2,
                              double[:,::1] alpha, double auto_z, double[::1] autoprior, double[::1] pk, double[::1] kcdf,
                              double[::1] prior_product, double[:,::1] B, double[:,::1] Bi, double[::1] u, long[::1] k, long pi,
                              double[::1] this_mu, double[::1] this_theta, double v):
    """
    Computes the Ordinary Least Square method for each model given by the design
    matrices X for the original data and estimates the probability of each model
    
    Parameters
    ----------
    data : double** (npoints, 3)
        Original data (x, y, std)   
    fg : double* (npoints)
        Value of the global function in each point
    fl : double* (npoints)
        Value of the local functions accumulated in each points
    xi : long
        Starting index of the partition
    xj : long
        Ending index of the partition
    kmax : long
        Maximum number of models
    nlocal_paramaters : long* (kmax)
        Number of local parameters for each model
    mu : double** (kmax, nlmax)
        *This is an output, original input will be overwritten* Mean values of the local parameters of each model
        obtained using the OLS method
    sigma : double** (kmax, nlmax)
        *This is an output, original input will be overwritten* Covariance matrix of the local parameters of each model
        obtained using the OLS method
    detS : double* (kmax)
        *This is an output, original input will be overwritten* Determinants of the covariance matrices
    epsilon : double* (kmax)
        *This is an output, original input will be overwritten* Residuals of each model
    X : double*** (kmax, npoints, nlmax)
        Design matrices of each model
    N : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Normal matrix (X^T*X)
    Q : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Inverse of normal matrix (N^-1)
    L : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Lower triangular matrix
    U : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Upper triangular matrix
    P : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Pivot matrix
    aux1 : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Auxiliar matrix for internal calculations
    aux2 : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Auxiliar matrix for internal calculations
    alpha : double** (kmax, kmax)
        *Temporary array, original input will be overwritten* Probability decomposition of each model
    auto_z : double
        Predefined prior for OLS estimation for all models
    autoprior : double* (kmax)
        *This is an output, original input will be overwritten* Autoprior of each model after OLS
    pk : double* (kmax)
        *This is an output, original input will be overwritten* Posterior probability of each model
    kcdf : double* (kmax)
        *This is an output, original input will be overwritten* Cumulative distribution of the models
        
    Returns
    -------
    valid : long
        If valid is equal to 1, the OLS and MCMC were performed successfully
    """
    cdef long x
    cdef long ki
    cdef long kj
    cdef long li
    cdef long lj
    cdef long lk
    cdef long ll
    cdef long nl
    cdef long nli
    cdef long nlj

    cdef double ep
    cdef double s2
    cdef double count
    
    cdef double detB
    cdef double sumpdf
    cdef double pdf
    cdef double curve_prob
    cdef double ppratio
    
    for ki in range(kmax):
        nl = nlocal_parameters[ki]
        for li in range(nl):
            for lj in range(nl):
                N[li][lj] = 0.
                for x in range(xi, xj+1):
                    N[li][lj] += X[ki][x][li]*X[ki][x][lj]
                    
        mat_zero(aux1, nl)
        mat_zero(aux2, nl)
        
        count = mat_LU(N, L, U, P, nl, aux1)
        mat_inverse_triangular_lower(L, aux1, nl)
        mat_inverse_triangular_upper(U, aux2, nl)
        
        mat_zero(Q, nl)
        
        for li in range(nl):
            for lj in range(li+1):
                Q[li][lj] = 0.
                for lk in range(nl):
                    for ll in range(nl):
                        Q[li][lj] += aux2[li][lk]*aux1[lk][ll]*P[ll][lj]
                Q[lj][li] = Q[li][lj]
        
        for li in range(nl):
            mu[ki][li] = 0.
            for lj in range(nl):
                for x in range(xi, xj+1):
                    mu[ki][li] += Q[li][lj]*X[ki][x][lj]*(data[x][1] - fg[x] - fl[x])
        
        epsilon[ki] = 0.
        s2 = 0.
        
        for x in range(xi, xj+1):
            ep = (data[x][1] - fg[x] - fl[x])
            for li in range(nl):
                ep -=  X[ki][x][li]*mu[ki][li]
            s2 += ep*ep
            epsilon[ki] += (ep*ep)/(2.*data[x][2]*data[x][2])
        
        s2 = s2/(<double> (xj - xi + 1 - nl))
        
        autoprior[ki] = 1.
    
        mat_zero(S[ki], nl)
        for li in range(nl):
            for lj in range(li+1):
                S[ki][li][lj] = s2*Q[li][lj]
                S[ki][lj][li] = S[ki][li][lj]
            sigma[ki][li] = sqrt(S[ki][li][li])
            autoprior[ki] *= 2.*auto_z*sigma[ki][li]
        
        detS[ki] = mat_LU(S[ki], L, U, P, nl, aux1)
    
    mat_zero(alpha, kmax)
    for ki in range(kmax):
        nli = nlocal_parameters[ki]
        alpha[ki][ki] = 1.
        for kj in range(ki+1, kmax):
            nlj = nlocal_parameters[kj]
            
            if detS[kj] > 0.:
                alpha[ki][kj] = exp(epsilon[kj] - epsilon[ki])*sqrt(pow(2.0 * M_PI, nli - nlj) * detS[ki]/detS[kj])
                alpha[ki][kj] *= autoprior[kj]/autoprior[ki]
            else:
                alpha[ki][kj] = INFINITY
            
            if alpha[ki][kj] == 0.:
                alpha[kj][ki] = INFINITY
            else:
                alpha[kj][ki] = 1./alpha[ki][kj]
    
    for ki in range(kmax):
        if ki == 0:
            kcdf[ki] = 0.
        else:
            kcdf[ki] = kcdf[ki-1]
            
        pk[ki] = 0.
        for kj in range(kmax):
            pk[ki] += alpha[kj][ki]
            
        pk[ki] = 1./pk[ki]
        kcdf[ki] += pk[ki]
        
    # Now update the prior product for each order
    for ki in range(kmax):
        prior_product[ki] = 1.
        nl = nlocal_parameters[ki]
        for li in range(nl):
            prior_product[ki] *= 2.*auto_z*sigma[ki][li]
            
    # Resample partition
    # Choose the new sub-model
    ki = choose_interval(kcdf, kmax, v)
    if ki < 0:
        return -1.
    k[pi] = ki
    
    # Compute the Cholesky decomposition of the given submodel
    nl = nlocal_parameters[ki]
    status = Cholesky_Decomposition(S[ki], B, nl)
    if status < 0:
        return -1.
    
    # Compute inverse Cholesky decomposition
    mat_inverse_triangular_lower(B, Bi, nl)
    
    # Save mean values and compute the inverse of the Cholesky decomposition matrix
    detB = 1.
    for li in range(nl):
        detB *= B[li][li]
        this_mu[li] = mu[ki][li]
        
    # Sample a random curve using a multivariate normal distribution from the best fit curve
    for li in range(nl):
        this_theta[li] = this_mu[li]
        for lj in range(nl):
            this_theta[li] += B[li][lj]*u[lj]
            
    # Evaluate the multivariate normal distribution pdf
    sumpdf = 0.
    for li in range(nl):
        pdf = 0.
        for lj in range(nl):
            pdf += Bi[li][lj]*(this_theta[lj] - this_mu[lj])
        sumpdf += pdf*pdf
        
    curve_prob = 1./(pow(2.0*M_PI, <double>(nl)/2.0) * detB) * exp(-0.5*sumpdf)
    
    ppratio = pk[ki]/(prior_product[ki]*curve_prob)
    
    return ppratio

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double compute_fixed_init(double[:,::1] data, double[::1] fg, double[::1] fl, long xi, long xj, long kmax,
                              long[::1] nlocal_parameters, double[:,::1] mu, double[:,::1] sigma, double[:,:,::1] S,
                              double[::1] detS, double[::1] epsilon, double[:,:,::1] X, double[:,::1] N, double[:,::1] Q,
                              double[:,::1] L, double[:,::1] U, double[:,::1] P, double[:,::1] aux1, double[:,::1] aux2,
                              double[:,::1] alpha, double auto_z, double[::1] autoprior, double[::1] pk, double[::1] kcdf,
                              double[::1] prior_product, double[:,::1] B, double[:,::1] Bi, double[::1] u, long[::1] k, long pi,
                              double[::1] this_mu, double[::1] this_theta, long this_k, long use_theta):
    """
    Computes the Ordinary Least Square method for each model given by the design
    matrices X for the original data and estimates the probability of each model
    
    Parameters
    ----------
    data : double** (npoints, 3)
        Original data (x, y, std)   
    fg : double* (npoints)
        Value of the global function in each point
    fl : double* (npoints)
        Value of the local functions accumulated in each points
    xi : long
        Starting index of the partition
    xj : long
        Ending index of the partition
    kmax : long
        Maximum number of models
    nlocal_paramaters : long* (kmax)
        Number of local parameters for each model
    mu : double** (kmax, nlmax)
        *This is an output, original input will be overwritten* Mean values of the local parameters of each model
        obtained using the OLS method
    sigma : double** (kmax, nlmax)
        *This is an output, original input will be overwritten* Covariance matrix of the local parameters of each model
        obtained using the OLS method
    detS : double* (kmax)
        *This is an output, original input will be overwritten* Determinants of the covariance matrices
    epsilon : double* (kmax)
        *This is an output, original input will be overwritten* Residuals of each model
    X : double*** (kmax, npoints, nlmax)
        Design matrices of each model
    N : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Normal matrix (X^T*X)
    Q : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Inverse of normal matrix (N^-1)
    L : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Lower triangular matrix
    U : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Upper triangular matrix
    P : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Pivot matrix
    aux1 : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Auxiliar matrix for internal calculations
    aux2 : double** (nlmax, nlmax)
        *Temporary array, original input will be overwritten* Auxiliar matrix for internal calculations
    alpha : double** (kmax, kmax)
        *Temporary array, original input will be overwritten* Probability decomposition of each model
    auto_z : double
        Predefined prior for OLS estimation for all models
    autoprior : double* (kmax)
        *This is an output, original input will be overwritten* Autoprior of each model after OLS
    pk : double* (kmax)
        *This is an output, original input will be overwritten* Posterior probability of each model
    kcdf : double* (kmax)
        *This is an output, original input will be overwritten* Cumulative distribution of the models
        
    Returns
    -------
    valid : long
        If valid is equal to 1, the OLS and MCMC were performed successfully
    """
    cdef long x
    cdef long ki
    cdef long kj
    cdef long li
    cdef long lj
    cdef long lk
    cdef long ll
    cdef long nl
    cdef long nli
    cdef long nlj

    cdef double ep
    cdef double s2
    cdef double count
    
    cdef double detB
    cdef double sumpdf
    cdef double pdf
    cdef double curve_prob
    cdef double ppratio
    
    for ki in range(kmax):
        nl = nlocal_parameters[ki]
        for li in range(nl):
            for lj in range(nl):
                N[li][lj] = 0.
                for x in range(xi, xj+1):
                    N[li][lj] += X[ki][x][li]*X[ki][x][lj]
                    
        mat_zero(aux1, nl)
        mat_zero(aux2, nl)
        
        count = mat_LU(N, L, U, P, nl, aux1)
        mat_inverse_triangular_lower(L, aux1, nl)
        mat_inverse_triangular_upper(U, aux2, nl)
        
        mat_zero(Q, nl)
        
        for li in range(nl):
            for lj in range(li+1):
                Q[li][lj] = 0.
                for lk in range(nl):
                    for ll in range(nl):
                        Q[li][lj] += aux2[li][lk]*aux1[lk][ll]*P[ll][lj]
                Q[lj][li] = Q[li][lj]
        
        for li in range(nl):
            mu[ki][li] = 0.
            for lj in range(nl):
                for x in range(xi, xj+1):
                    mu[ki][li] += Q[li][lj]*X[ki][x][lj]*(data[x][1] - fg[x] - fl[x])
        
        epsilon[ki] = 0.
        s2 = 0.
        
        for x in range(xi, xj+1):
            ep = (data[x][1] - fg[x] - fl[x])
            for li in range(nl):
                ep -=  X[ki][x][li]*mu[ki][li]
            s2 += ep*ep
            epsilon[ki] += (ep*ep)/(2.*data[x][2]*data[x][2])
        
        s2 = s2/(<double> (xj - xi + 1 - nl))
        
        autoprior[ki] = 1.
    
        mat_zero(S[ki], nl)
        for li in range(nl):
            for lj in range(li+1):
                S[ki][li][lj] = s2*Q[li][lj]
                S[ki][lj][li] = S[ki][li][lj]
            sigma[ki][li] = sqrt(S[ki][li][li])
            autoprior[ki] *= 2.*auto_z*sigma[ki][li]
        
        detS[ki] = mat_LU(S[ki], L, U, P, nl, aux1)
    
    mat_zero(alpha, kmax)
    for ki in range(kmax):
        nli = nlocal_parameters[ki]
        alpha[ki][ki] = 1.
        for kj in range(ki+1, kmax):
            nlj = nlocal_parameters[kj]
            
            if detS[kj] > 0.:
                alpha[ki][kj] = exp(epsilon[kj] - epsilon[ki])*sqrt(pow(2.0 * M_PI, nli - nlj) * detS[ki]/detS[kj])
                alpha[ki][kj] *= autoprior[kj]/autoprior[ki]
            else:
                alpha[ki][kj] = INFINITY
            
            if alpha[ki][kj] == 0.:
                alpha[kj][ki] = INFINITY
            else:
                alpha[kj][ki] = 1./alpha[ki][kj]
    
    for ki in range(kmax):
        if ki == 0:
            kcdf[ki] = 0.
        else:
            kcdf[ki] = kcdf[ki-1]
            
        pk[ki] = 0.
        for kj in range(kmax):
            pk[ki] += alpha[kj][ki]
            
        pk[ki] = 1./pk[ki]
        kcdf[ki] += pk[ki]
        
    # Now update the prior product for each order
    for ki in range(kmax):
        prior_product[ki] = 1.
        nl = nlocal_parameters[ki]
        for li in range(nl):
            prior_product[ki] *= 2.*auto_z*sigma[ki][li]
            
    # Resample partition
    # Choose the new sub-model
    ki = this_k
    k[pi] = ki
    
    # Compute the Cholesky decomposition of the given submodel
    nl = nlocal_parameters[ki]
    status = Cholesky_Decomposition(S[ki], B, nl)
    if status < 0:
        return -1.
    
    # Compute inverse Cholesky decomposition
    mat_inverse_triangular_lower(B, Bi, nl)
    
    # Save mean values and compute the inverse of the Cholesky decomposition matrix
    detB = 1.
    for li in range(nl):
        detB *= B[li][li]
        this_mu[li] = mu[ki][li]
        
    # Sample a random curve using a multivariate normal distribution from the best fit curve
    if not use_theta:
        for li in range(nl):
            this_theta[li] = this_mu[li]
            for lj in range(nl):
                this_theta[li] += B[li][lj]*u[lj]
            
    # Evaluate the multivariate normal distribution pdf
    sumpdf = 0.
    for li in range(nl):
        pdf = 0.
        for lj in range(nl):
            pdf += Bi[li][lj]*(this_theta[lj] - this_mu[lj])
        sumpdf += pdf*pdf
        
    curve_prob = 1./(pow(2.0*M_PI, <double>(nl)/2.0) * detB) * exp(-0.5*sumpdf)
    
    ppratio = pk[ki]/(prior_product[ki]*curve_prob)
    
    return ppratio

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double update_partition_single_model(double[:,::1] data, double[::1] fg, double[::1] fl, long xi, long xj,
                                           long nl, double[::1] mu, double[::1] sigma,
                                           double[:,::1] S, double[:,::1] X, double[:,::1] N, double[:,::1] Q,
                                           double[:,::1] L, double[:,::1] U, double[:,::1] P, double[:,::1] aux1,
                                           double[:,::1] aux2, double auto_z,
                                           double[::1] prior_product, double[:,::1] B, double[:,::1] Bi,
                                           double[::1] u, double[::1] this_mu, double[::1] this_theta):
    
    cdef long x
    cdef long li
    cdef long lj
    cdef long lk
    cdef long ll

    cdef double epsilon
    cdef double ep
    cdef double s2
    cdef double count
    
    cdef double detS
    cdef double detB
    cdef double sumpdf
    cdef double pdf
    cdef double curve_prob
    cdef double ppratio
    
    
    for li in range(nl):
        for lj in range(nl):
            N[li][lj] = 0.
            for x in range(xi, xj):
                N[li][lj] += X[x][li]*X[x][lj]
                
    mat_zero(aux1, nl)
    mat_zero(aux2, nl)
    
    count = mat_LU(N, L, U, P, nl, aux1)
    mat_inverse_triangular_lower(L, aux1, nl)
    mat_inverse_triangular_upper(U, aux2, nl)
    
    mat_zero(Q, nl)
    
    for li in range(nl):
        for lj in range(li+1):
            Q[li][lj] = 0.
            for lk in range(nl):
                for ll in range(nl):
                    Q[li][lj] += aux2[li][lk]*aux1[lk][ll]*P[ll][lj]
            Q[lj][li] = Q[li][lj]
    
    for li in range(nl):
        mu[li] = 0.
        for lj in range(nl):
            for x in range(xi, xj):
                mu[li] += Q[li][lj]*X[x][lj]*(data[x][1] - fg[x] - fl[x])
    
    epsilon = 0.
    s2 = 0.
    
    for x in range(xi, xj):
        ep = (data[x][1] - fg[x] - fl[x])
        for li in range(nl):
            ep -=  X[x][li]*mu[li]
        s2 += ep*ep
        epsilon += (ep*ep)/(2.*data[x][2]*data[x][2])
    
    s2 = s2/(<double> (xj - xi - nl))

    mat_zero(S, nl)
    for li in range(nl):
        for lj in range(li+1):
            S[li][lj] = s2*Q[li][lj]
            S[lj][li] = S[li][lj]
        sigma[li] = sqrt(S[li][li])
    
    detS = mat_LU(S, L, U, P, nl, aux1)
    
    alpha = 1.
        
    # Now update the prior product
    prior_product[0] = 1.
    for li in range(nl):
        prior_product[0] *= 2.*auto_z*sigma[li]
            
    # Compute the Cholesky decomposition
    status = Cholesky_Decomposition(S, B, nl)
    if status < 0:
        return -1.
    
    # Compute inverse Cholesky decomposition
    mat_inverse_triangular_lower(B, Bi, nl)
    
    # Save mean values and compute the inverse of the Cholesky decomposition matrix
    detB = 1.
    for li in range(nl):
        detB *= B[li][li]
        this_mu[li] = mu[li]
        
    # Sample a random curve using a multivariate normal distribution from the best fit curve
    for li in range(nl):
        this_theta[li] = this_mu[li]
        for lj in range(nl):
            this_theta[li] += B[li][lj]*u[lj]
            
    # Evaluate the multivariate normal distribution pdf
    sumpdf = 0.
    for li in range(nl):
        pdf = 0.
        for lj in range(nl):
            pdf += Bi[li][lj]*(this_theta[lj] - this_mu[lj])
        sumpdf += pdf*pdf
        
    curve_prob = 1./(pow(2.0*M_PI, <double>(nl)/2.0) * detB) * exp(-0.5*sumpdf)
    
    ppratio = 1./(prior_product[0]*curve_prob)
    
    return ppratio
