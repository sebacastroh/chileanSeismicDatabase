# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 19:21:07 2021

@author: sebac
"""

import math
import numpy as np

# Single partition models
from models.single_partition.single_order_independent_uniform import sp_so_iu
from models.single_partition.single_order_random_walk_gauss import sp_so_rwg
from models.single_partition.single_order_gauss_ols import sp_so_go

# # from models.single_partition.multiple_order_independent_uniform import sp_mo_iu
# # from models.single_partition.multiple_order_random_walk_gauss import sp_mo_rwg

from models.single_partition.multiple_order_gauss_ols import sp_mo_go

# Multiple partition models
from models.multiple_partition.single_order_independent_uniform import mp_so_iu
from models.multiple_partition.single_order_random_walk_gauss import mp_so_rwg
from models.multiple_partition.single_order_gauss_ols_pulse import mp_so_gop
from models.multiple_partition.single_order_gauss_ols_heaviside import mp_so_goh

# # from models.multiple_partition.multiple_order_independent_uniform import mp_mo_iu
# # from models.multiple_partition.multiple_order_random_walk_gauss import mp_mo_rwg

from models.multiple_partition.multiple_order_gauss_ols_pulse import mp_mo_gop
from models.multiple_partition.multiple_order_gauss_ols_heaviside import mp_mo_goh

from models.multiple_partition.single_order_connected_lines import mp_so_connected_lines

from results import result1D

def _engine(model, datasets, burnin, total, min_partitions, max_partitions,
            kmax, samples, pd, options, weights, nglobal_parameters,
            global_parameters, nlocal_parameters, local_parameters,
            computeDesignMatrices, loglikelihood_function, global_function,
            design_matrices, local_function, sampling_function, nprocesses,
            fixed_init, seed):
    
    nprs = np.random.default_rng(seed)
        
    ndatasets = len(datasets)
    npoints = []
    for dataset in datasets:
        npoints.append(len(dataset))
    npoints = np.array(npoints)
    
    xmin = np.inf
    xmax = -np.inf
    
    for di in range(ndatasets):
        xmin = min(xmin, datasets[di][:,0].min())
        xmax = max(xmax, datasets[di][:,0].max())
        
    # Create results object
    results = result1D(ndatasets, kmax, max_partitions, burnin, total,
                       nglobal_parameters, nlocal_parameters,
                       samples, xmin, xmax)
    
    # Create model objects
    current = model(npoints, ndatasets, xmin, xmax, min_partitions, max_partitions,
                 kmax, samples, pd, options, weights, nglobal_parameters,
                 global_parameters, nlocal_parameters, local_parameters,
                 computeDesignMatrices, loglikelihood_function, global_function,
                 design_matrices, local_function, sampling_function, fixed_init, nprs)
    
    proposed = model(npoints, ndatasets, xmin, xmax, min_partitions, max_partitions,
                 kmax, samples, pd, options, weights, nglobal_parameters,
                 global_parameters, nlocal_parameters, local_parameters,
                 computeDesignMatrices, loglikelihood_function, global_function,
                 design_matrices, local_function, sampling_function, fixed_init, nprs)
    
    # Initialize model
    if fixed_init is None or len(fixed_init) == 0:
        current.initialize(datasets)
    else:
        success = current.initialize_fixed(datasets)
        if not success:
            raise 'Fixed init model failed'
            
    
    if sampling_function is not None:
        x = np.linspace(xmin, xmax, samples)
        y = []
        for di in range(ndatasets):
            y.append(current.sample(x, di))
    else:
        y = None
        
    # Compute misfit
    current.misfit(datasets)
    
    for i in range(total):
        # Select random process
        if options is None:
            process = 0
        else:
            process = nprs.choice(options, p=weights)
        
        # Clone model
        proposed.clone(current)
        
        # Perturb proposed model
        success, process_prob = proposed.perturb(process, datasets)
        if not success:
            results._sample(current, y, process, False)
            continue
        
        # Compute misfit
        proposed.misfit(datasets)
        
        # Compute acceptance
        # print('---')
        # print(process)
        # print(process_prob)
        prob = math.log(process_prob) + current.likelihood - proposed.likelihood
        
        # Sample random value, if u < alpha then the proposed model is accepted
        # and then the current model is updated
        u = math.log(nprs.random())
        if u < prob:
            if sampling_function is not None:
                for di in range(ndatasets):
                    y[di] = proposed.sample(x, di)
            current.clone(proposed)
            results._sample(current, y, process, True)
        else:
            results._sample(current, y, process, False)
    
    # Assemble results
    results._assemble(nprocesses)
    
    return results

def independent_metropolis_hastings_uniform(dataset, burnin, total, samples,
                                           nparameters, parameters, loglikelihood_function,
                                           sampling_function, fixed_init=None, seed=None):
    
    min_partitions = 1
    max_partitions = 1
    kmax = 1
    pd = 0.
    options = None
    weights = None
    nglobal_parameters = 0
    global_parameters = None
    nlocal_parameters = nparameters
    local_parameters = parameters
    computeDesignMatrices = False
    global_function = None
    design_matrices = None
    local_function = None
    nprocesses = 1
    
    model = sp_so_iu
    results = _engine(model, [dataset], burnin, total, min_partitions, max_partitions,
                      kmax, samples, pd, options, weights, nglobal_parameters,
                      global_parameters, nlocal_parameters, local_parameters,
                      computeDesignMatrices, loglikelihood_function, global_function,
                      design_matrices, local_function, sampling_function, nprocesses,
                      fixed_init, seed)
    
    return results

def metropolis_hastings_random_walk_gauss(dataset, burnin, total, samples,
                                         nparameters, parameters, loglikelihood_function,
                                         sampling_function, fixed_init=None, seed=None):
    
    min_partitions = 1
    max_partitions = 1
    kmax = 1
    pd = 0.
    options = None
    weights = None
    nglobal_parameters = 0
    global_parameters = None
    nlocal_parameters = nparameters
    local_parameters = parameters
    computeDesignMatrices = False
    global_function = None
    design_matrices = None
    local_function = None
    nprocesses = 1
    
    model = sp_so_rwg
    results = _engine(model, [dataset], burnin, total, min_partitions, max_partitions,
                      kmax, samples, pd, options, weights, nglobal_parameters,
                      global_parameters, nlocal_parameters, local_parameters,
                      computeDesignMatrices, loglikelihood_function, global_function,
                      design_matrices, local_function, sampling_function, nprocesses,
                      fixed_init, seed)
    
    return results

def metropolis_hastings_ols(dataset, burnin, total, samples,
                            nlocal_parameters,
                            design_matrix, loglikelihood_function,
                            sampling_function=None,
                            computeDesignMatrix=False,
                            nglobal_parameters=0, global_parameters=None,
                            global_function=None, weights=None,
                            fixed_init=None, seed=None):
    
    min_partitions = 1
    max_partitions = 1
    kmax = 1
    pd = 0.
    local_parameters = None
    
    if nglobal_parameters == 0:
        options = None
        nprocesses = 1
    else:
        nprocesses = 2
        options = np.arange(2)
        if weights is None:
            weights = 0.5*np.ones(2)
    
    local_function = None
    
    model = sp_so_go
    results = _engine(model, [dataset], burnin, total, min_partitions, max_partitions,
                      kmax, samples, pd, options, weights, nglobal_parameters,
                      global_parameters, nlocal_parameters, local_parameters,
                      computeDesignMatrix, loglikelihood_function, global_function,
                      design_matrix, local_function, sampling_function, nprocesses,
                      fixed_init, seed)
    
    return results

def single_partition_multiple_order_ols(dataset, burnin, total, samples, nlocal_parameters,
                                        design_matrices, loglikelihood_function,
                                        sampling_function=None,
                                        computeDesignMatrices=False,
                                        nglobal_parameters=0, global_parameters=None,
                                        global_function=None, weights=None,
                                        fixed_init=None, seed=None):
    
    min_partitions = 1
    max_partitions = 1
    kmax = len(nlocal_parameters)
    pd = 0.
    
    if nglobal_parameters == 0:
        options = None
        nprocesses = 1
    else:
        nprocesses = 2
        options = np.arange(2)
        if weights is None:
            weights = 0.5*np.ones(2)
    
    local_parameters = None
    local_function = None
    
    model = sp_mo_go
    results = _engine(model, [dataset], burnin, total, min_partitions, max_partitions,
                      kmax, samples, pd, options, weights, nglobal_parameters,
                      global_parameters, nlocal_parameters, local_parameters,
                      computeDesignMatrices, loglikelihood_function, global_function,
                      design_matrices, local_function, sampling_function, nprocesses,
                      fixed_init, seed)
    
    return results

def single_partition_polynomial_fit(dataset, burnin, total, samples, kmax, seed=None):
    
    if kmax == 0:
        nlocal_parameters = 1
    else:
        nlocal_parameters = np.arange(kmax + 1) + 1
        
    nlmax = kmax + 1
    
    nglobal_parameters = 0
    global_parameters = None
    global_function = None
    weights = None
    fixed_init = None
    
    npoints = len(dataset)
    design_matrices = np.ones((nlmax, npoints, nlmax))
    
    for ki in range(kmax+1):
        for li in range(ki+1):
            design_matrices[ki,:,li] = dataset[:,0]**li
    
    computeDesignMatrices = False
    
    def sampling_function(x, ki, local_values):
        X = np.zeros((len(x), len(local_values)))
        for li in range(ki+1):
            X[:,li] = x**li
        y = np.dot(X, local_values)
        return y
    
    def loglikelihood_function(ki, local_values, x, data, sigma):
        y = sampling_function(x, ki, local_values)
        phi = np.sum(((y - data)/(2.*sigma))**2)
        
        return phi
    
    results = single_partition_multiple_order_ols(dataset, burnin, total, samples, nlocal_parameters,
                                                  design_matrices, loglikelihood_function,
                                                  sampling_function, computeDesignMatrices,
                                                  nglobal_parameters, global_parameters,
                                                  global_function, weights,
                                                  fixed_init, seed)
    
    return results

def multiple_partition_single_order_independent_uniform(datasets, burnin, total, min_partitions,
                                                        max_partitions, samples, pd,
                                                        nlocal_parameters, local_parameters,
                                                        loglikelihood_function, sampling_function=None,
                                                        nglobal_parameters=0, global_parameters=None,
                                                        weights=None, fixed_init=None, seed=None):
    
    kmax = 1
    local_function = None
    
    if nglobal_parameters == 0:
        nprocesses = 5
    else:
        nprocesses = 6
    
    options = np.arange(nprocesses)
    
    if weights is None:
        weights = np.ones(nprocesses)/nprocesses
        
    design_matrix = None
    computeDesignMatrix = False
    global_function = None
    
    model = mp_so_iu
    
    results = _engine(model, datasets, burnin, total, min_partitions, max_partitions,
            kmax, samples, pd, options, weights, nglobal_parameters,
            global_parameters, nlocal_parameters, local_parameters,
            computeDesignMatrix, loglikelihood_function, global_function,
            design_matrix, local_function, sampling_function, nprocesses,
            fixed_init, seed)
    
    return results

def multiple_partition_single_order_random_walk_gaussian(datasets, burnin, total, min_partitions,
                                                         max_partitions, samples, pd,
                                                         nlocal_parameters, local_parameters,
                                                         loglikelihood_function, sampling_function=None,
                                                         nglobal_parameters=0, global_parameters=None,
                                                         weights=None, fixed_init=None, seed=None):
    
    kmax = 1
    local_function = None
    
    if nglobal_parameters == 0:
        nprocesses = 5
    else:
        nprocesses = 6
    
    options = np.arange(nprocesses)
    
    if weights is None:
        weights = np.ones(nprocesses)/nprocesses
        
    design_matrix = None
    computeDesignMatrix = False
    global_function = None
    
    model = mp_so_rwg
    
    results = _engine(model, datasets, burnin, total, min_partitions, max_partitions,
            kmax, samples, pd, options, weights, nglobal_parameters,
            global_parameters, nlocal_parameters, local_parameters,
            computeDesignMatrix, loglikelihood_function, global_function,
            design_matrix, local_function, sampling_function, nprocesses,
            fixed_init, seed)
    
    return results

def multiple_partition_single_order_pulse_ols(datasets, burnin, total, min_partitions,
                                              max_partitions, samples, pd,
                                              nlocal_parameters, design_matrix,
                                              loglikelihood_function, sampling_function=None,
                                              nglobal_parameters=0, global_parameters=None,
                                              global_function=None,
                                              weights=None, fixed_init=None, seed=None):
    
    kmax = 1
    local_parameters = None
    local_function = None
    
    if nglobal_parameters == 0:
        nprocesses = 5
    else:
        nprocesses = 6
    
    options = np.arange(nprocesses)
    
    if weights is None:
        weights = np.ones(nprocesses)/nprocesses
        
    if callable(design_matrix):
        computeDesignMatrix = True
    else:
        computeDesignMatrix = False
    
    
    model = mp_so_gop
    
    results = _engine(model, datasets, burnin, total, min_partitions, max_partitions,
            kmax, samples, pd, options, weights, nglobal_parameters,
            global_parameters, nlocal_parameters, local_parameters,
            computeDesignMatrix, loglikelihood_function, global_function,
            design_matrix, local_function, sampling_function, nprocesses,
            fixed_init, seed)
    
    return results

def multiple_partition_single_order_heaviside_ols(datasets, burnin, total, min_partitions,
                                                  max_partitions, samples, pd,
                                                  nlocal_parameters, design_matrix,
                                                  local_function, loglikelihood_function,
                                                  sampling_function=None, nglobal_parameters=0,
                                                  global_parameters=None, global_function=None,
                                                  weights=None, fixed_init=None, seed=None):
    
    kmax = 1
    local_parameters = None
    
    if nglobal_parameters == 0:
        nprocesses = 5
    else:
        nprocesses = 6
    
    options = np.arange(nprocesses)
    
    if weights is None:
        weights = np.ones(nprocesses)/nprocesses
        
    if callable(design_matrix):
        computeDesignMatrix = True
    else:
        computeDesignMatrix = False
    
    
    model = mp_so_goh
    
    results = _engine(model, datasets, burnin, total, min_partitions, max_partitions,
            kmax, samples, pd, options, weights, nglobal_parameters,
            global_parameters, nlocal_parameters, local_parameters,
            computeDesignMatrix, loglikelihood_function, global_function,
            design_matrix, local_function, sampling_function, nprocesses,
            fixed_init, seed)
    
    return results

def multiple_partition_multiple_order_pulse_ols(datasets, burnin, total, min_partitions,
                                                max_partitions, samples, pd, kmax,
                                                nlocal_parameters, design_matrices,
                                                loglikelihood_function, sampling_function=None,
                                                nglobal_parameters=0, global_parameters=None,
                                                global_function=None,
                                                weights=None, fixed_init=None, seed=None):
    
    local_parameters = None
    local_function = None
    
    if nglobal_parameters == 0:
        nprocesses = 5
    else:
        nprocesses = 6
    
    options = np.arange(nprocesses)
    
    if weights is None:
        weights = np.ones(nprocesses)/nprocesses
        
    if callable(design_matrices):
        computeDesignMatrices = True
    else:
        computeDesignMatrices = False
    
    
    model = mp_mo_gop
    
    results = _engine(model, datasets, burnin, total, min_partitions, max_partitions,
            kmax, samples, pd, options, weights, nglobal_parameters,
            global_parameters, nlocal_parameters, local_parameters,
            computeDesignMatrices, loglikelihood_function, global_function,
            design_matrices, local_function, sampling_function, nprocesses,
            fixed_init, seed)
    
    return results

def multiple_partition_polynomial_fit(datasets, burnin, total, min_partitions,
                                      max_partitions, samples, pd, kmax, seed=None):
    
    nlmax = kmax + 1
    nlocal_parameters = np.arange(nlmax) + 1
    
    design_matrices = []
    for di, dataset in enumerate(datasets):
        npoints = len(dataset)
        X = np.ones((nlmax, npoints, nlmax))
        
        for ki in range(kmax+1):
            for li in range(ki+1):
                X[ki,:,li] = dataset[:,0]**li
                
        design_matrices.append(X)
    
    def sampling_function(x, boundaries, kvalues, local_values, nlocal_parameters):
        
        y = np.zeros_like(x)
        for i in range(len(boundaries)-1):
            pos = np.where((x >= boundaries[i]) & (x <= boundaries[i+1]))
            y[pos] = np.polyval(local_values[i,:nlocal_parameters[ki]][::-1], x[pos])
        
        return y
    
    def loglikelihood_function(boundaries, kvalues, local_values, nlocal_parameters, x, data, sigma):
        
        y = sampling_function(x, boundaries, kvalues, local_values, nlocal_parameters)
        phi = np.sum(((y - data)/(2.*sigma))**2)
        
        return phi
    
    nglobal_parameters = 0
    global_parameters = None
    global_function = None
    weights = None
    fixed_init = None
    
    results = multiple_partition_multiple_order_pulse_ols(datasets, burnin, total, min_partitions,
                                                          max_partitions, samples, pd, kmax,
                                                          nlocal_parameters, design_matrices,
                                                          loglikelihood_function, sampling_function,
                                                          nglobal_parameters, global_parameters,
                                                          global_function,
                                                          weights, fixed_init, seed)
    
    return results

def multiple_partition_multiple_order_heaviside_ols(datasets, burnin, total, min_partitions,
                                                    max_partitions, samples, pd, kmax,
                                                    nlocal_parameters, design_matrices, local_function,
                                                    loglikelihood_function, sampling_function=None,
                                                    nglobal_parameters=0, global_parameters=None,
                                                    global_function=None,
                                                    weights=None, fixed_init=None, seed=None):
    
    local_parameters = None
    
    if nglobal_parameters == 0:
        nprocesses = 5
    else:
        nprocesses = 6
    
    options = np.arange(nprocesses)
    
    if weights is None:
        weights = np.ones(nprocesses)/nprocesses
        
    if callable(design_matrices):
        computeDesignMatrices = True
    else:
        computeDesignMatrices = False
    
    
    model = mp_mo_goh
    
    results = _engine(model, datasets, burnin, total, min_partitions, max_partitions,
            kmax, samples, pd, options, weights, nglobal_parameters,
            global_parameters, nlocal_parameters, local_parameters,
            computeDesignMatrices, loglikelihood_function, global_function,
            design_matrices, local_function, sampling_function, nprocesses,
            fixed_init, seed)
    
    return results

def multiple_partition_connected_lines(datasets, burnin, total, min_partitions,
                                       max_partitions, samples, pd, weights=None,
                                       seed=None):
    
    kmax = 1
    nlocal_parameters = 2
    local_parameters = None
    computeDesignMatrices = False
    design_matrices = None
    local_function = None
    
    nprocesses = 5
    options = np.arange(nprocesses)
    
    nglobal_parameters = 0
    global_parameters = None
    global_function = None
    if samples > 0:
        nsamples = samples
        sampling_function = True
    else:
        nsamples = 1
        sampling_function = None
    loglikelihood_function = None
    fixed_init = None
    
    model = mp_so_connected_lines
    
    results = _engine(model, datasets, burnin, total, min_partitions, max_partitions,
            kmax, nsamples, pd, options, weights, nglobal_parameters,
            global_parameters, nlocal_parameters, local_parameters,
            computeDesignMatrices, loglikelihood_function, global_function,
            design_matrices, local_function, sampling_function, nprocesses,
            fixed_init, seed)
    
    return results
