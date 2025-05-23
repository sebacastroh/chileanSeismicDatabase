# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 20:23:55 2021

@author: sebac
"""

import utils
import numpy as np

class sp_mo_go:
    
    def __init__(self, npoints, ndatasets, xmin, xmax, min_partitions, max_partitions,
                 kmax, samples, pd, options, weights, nglobal_parameters,
                 global_parameters, nlocal_parameters, local_parameters,
                 computeDesignMatrices, loglikelihood_function, global_function,
                 design_matrices, local_function, sampling_function, fixed_init, nprs):
        
        # Dataset parameters
        self.npoints = npoints[0]
        self.xmin = xmin
        self.xmax = xmax
        
        # Models parameters
        self.nglobal_parameters = nglobal_parameters
        self.kmax = kmax
        self.nlocal_parameters = nlocal_parameters.copy()
        self.nlmax = self.nlocal_parameters.max()
        
        if self.nglobal_parameters > 0:
            self.global_parameters = global_parameters.copy()
        
        # Model values
        self.theta_g = np.empty(self.nglobal_parameters)
        
        # Design matrices
        self.computeDesignMatrices = computeDesignMatrices
        self.X = np.zeros((self.kmax, self.npoints, self.nlmax))
        self.design_matrices = design_matrices
        
        # Current local models
        self.k = -np.ones(1, dtype=int) # Selected model
        self.theta = np.zeros(self.nlmax) # Sampled values of local parameters
        self.mu = np.zeros(self.nlmax) # Mean values of local parameters
        self.B = np.zeros((self.nlmax, self.nlmax)) # Cholesky decomposition of covariance matrices
        self.Bi = np.zeros((self.nlmax, self.nlmax)) # Inverses of the Cholesky decomposition of covariance matrices
        self.pk = np.zeros(self.kmax) # Prior probability of sub-models using OLS
        self.kcdf = np.zeros(self.kmax) # Cumulative prior probabilities
        self.prior_product = np.zeros(self.kmax)
        self.ppratio = 0
        
        # RJMCMC
        self.likelihood = 0.
        
        # Global function values
        self.gf_values = np.zeros(self.npoints)
        self.global_function = global_function
        
        # For internal computations
        self._u = np.zeros(self.nlmax) # Random vector
        self._mu = np.zeros((self.kmax, self.nlmax)) # Mean of local parameters
        self._sigma = np.ones((self.kmax, self.nlmax)) # Standard deviation of local parameters
        self._S = np.zeros((self.kmax, self.nlmax, self.nlmax)) # Covariance matrices of each local model
        self._detS = np.zeros(self.kmax) # Determinants of covariance matrices
        self._B = np.zeros((self.nlmax, self.nlmax)) # Cholesky decomposition of covariance matrix
        self._Bi = np.zeros((self.nlmax, self.nlmax)) # Inverse of the Cholesky decomposition of covariance matrix
        self._epsilon = np.zeros(self.kmax) # Residual of each local model
        self._N = np.zeros((self.nlmax, self.nlmax)) # Normal matrix
        self._Q = np.zeros((self.nlmax, self.nlmax)) # Cofactor matrix (N^T*N)^-1
        self._L = np.zeros((self.nlmax, self.nlmax)) # Lower triangular matrix
        self._U = np.zeros((self.nlmax, self.nlmax)) # Upper triangular matrix
        self._P = np.zeros((self.nlmax, self.nlmax)) # Permutation matrix
        self._aux1 = np.zeros((self.nlmax, self.nlmax)) # Auxiliar matrix
        self._aux2 = np.zeros((self.nlmax, self.nlmax)) # Auxiliar matrix
        self._autoprior = np.zeros(self.kmax)
        self._alpha = np.zeros((self.kmax, self.kmax)) # Acceptance probability matrix
        self._auto_z = 3.
        
        # Sampling function
        self.sampling_function = sampling_function
        
        # Likelihood function
        self.likelihood_function = loglikelihood_function
        
        self.nprs = nprs
    
    def initialize(self, datasets):
        
        dataset = datasets[0]
        # Sample global parameters if there is any
        if self.nglobal_parameters > 0:
            self.theta_g = self.nprs.random(self.nglobal_parameters)*\
                (self.global_parameters[:,1] - self.global_parameters[:,0]) + self.global_parameters[:,0]
                
            #Evaluate global function
            if callable(self.global_function):
                self.gf_values = self.global_function(dataset[:,0], self.theta_g)
        
        # Compute design matrices if corresponds
        if self.computeDesignMatrices:
            for ki in range(self.kmax):
                nl = self.nlocal_parameters[ki]
                self.X[ki,:,:nl] = self.design_matrices(dataset[:,0], self.theta_g, ki)
        else:
            self.X = self.design_matrices.copy()
        
        # Compute model for each segment
        self._u = self.nprs.normal(size=self.nlmax)
        self._update_partition(dataset)
        
    def initialize_fixed(self, datasets, fixed_init):
        
        pass
            
    def clone(self, other):
        
        # Model values
        self.theta_g = other.theta_g.copy()
        
        # Global function values
        self.gf_values = other.gf_values.copy()
        
        # Design matrices
        self.computeDesignMatrices = other.computeDesignMatrices
        self.X = other.X.copy()
        
        # Current local models
        self.theta = other.theta.copy()
        self.mu = other.mu.copy()
        self.B = other.B.copy()
        self.Bi = other.Bi.copy()
        self.pk = other.pk.copy()
        self.kcdf = other.kcdf.copy()
        self.prior_product = other.prior_product.copy()
        self.ppratio = other.ppratio
        
        # RJMCMC
        self.likelihood = other.likelihood
        
    
    def perturb(self, process, datasets):
        
        dataset = datasets[0]
        if process == 0:
            status, process_prob = self._propose_local_parameters(dataset)
        else:
            status, process_prob = self._propose_global_parameters(dataset)
        
        return status, process_prob
    
    def _propose_local_parameters(self, dataset):
        
        prob = 1./self.ppratio
        
        ki = self.k[0]
        nl = self.nlocal_parameters[ki]
        
        self._u[:nl] = np.dot(self.Bi[:nl,:nl], self.theta[:nl] - self.mu[:nl])
        
        if nl < self.nlmax:
            self._u[nl:] = self.nprs.normal(size=self.nlmax - nl)
                
        curve_prob = self._update_partition(dataset)
        
        if curve_prob <= 0:
            return 0, 0.
        
        prob *= curve_prob
        
        return 1, prob
    
    def _propose_global_parameters(self, dataset):
        
        if self.nglobal_parameters == 1:
            gi = 0
        else:
            gi = self.nprs.choice(self.nglobal_parameters)
            
        self.theta_g[gi] += self.nprs.normal(0., self.global_parameters[gi,2])
        
        if self.theta_g[gi] < self.global_parameters[gi,0] or \
            self.theta_g[gi] > self.global_parameters[gi,1]:
            return 0, 0.
        
        if callable(self.global_function):
            self.gf_values = self.global_function(dataset[:,0], self.theta_g)
        
        if self.computeDesignMatrices:
            for ki in range(self.kmax):
                nl = self.nlocal_parameters[ki]
                self.X[ki,:,:nl] = self.design_matrices(dataset[:,0], self.theta_g, ki)
        
        prob = 1./self.ppratio
        
        ki = self.k[0]
        nl = self.nlocal_parameters[ki]
        
        self._u[:nl] = np.dot(self.Bi[:nl,:nl], self.theta[:nl] - self.mu[:nl])
        
        if nl < self.nlmax:
            self._u[nl:] = self.nprs.normal(size=self.nlmax - nl)
        
        curve_prob = self._update_partition(dataset)
        if curve_prob <= 0:
            return 0, 0.
            
        prob *= curve_prob
        
        return 1, prob
    
    def misfit(self, datasets):
        
        dataset = datasets[0]
        
        if self.nglobal_parameters > 0:
            likelihood = self.likelihood_function(self.theta_g, self.k[0], self.theta, dataset[:,0], dataset[:,1], dataset[:,2])
        else:
            likelihood = self.likelihood_function(self.k[0], self.theta, dataset[:,0], dataset[:,1], dataset[:,2])
            
        self.likelihood = likelihood
        
    def sample(self, x, di):
        
        ki = self.k[0]
        nl = self.nlocal_parameters[ki]
        if self.nglobal_parameters > 0:
            return self.sampling_function(x, self.theta_g, ki, self.theta[:nl])
        else:
            return self.sampling_function(x, ki, self.theta[:nl])
    
    def _update_partition(self, dataset):
        # Compute OLS and determine model probabilities
        prob = utils.update_partition(dataset,
                                      self.gf_values,
                                      np.zeros(self.npoints),
                                      0,
                                      self.npoints-1,
                                      self.kmax,
                                      self.nlocal_parameters,
                                      self._mu,
                                      self._sigma,
                                      self._S,
                                      self._detS,
                                      self._epsilon,
                                      self.X,
                                      self._N,
                                      self._Q,
                                      self._L,
                                      self._U,
                                      self._P,
                                      self._aux1,
                                      self._aux2,
                                      self._alpha,
                                      self._auto_z,
                                      self._autoprior,
                                      self.pk,
                                      self.kcdf,
                                      self.prior_product,
                                      self.B,
                                      self.Bi,
                                      self._u,
                                      self.k,
                                      0,
                                      self.mu,
                                      self.theta,
                                      self.nprs.random())
        
        self.ppratio = prob
        
        return prob
