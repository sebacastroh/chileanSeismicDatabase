# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 21:21:36 2021

@author: sebac
"""
import utils
import numpy as np

class sp_so_go:
    
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
        self.nlocal_parameters = nlocal_parameters
        
        if self.nglobal_parameters > 0:
            self.global_parameters = global_parameters.copy()
        
        # Model values
        self.theta_g = np.empty(self.nglobal_parameters)
        
        # Design matrices
        self.computeDesignMatrices = computeDesignMatrices
        self.X = np.zeros((self.npoints, self.nlocal_parameters))
        self.design_matrices = design_matrices
        
        # Current local models
        self.theta = np.zeros(self.nlocal_parameters) # Sampled values of local parameters
        self.mu = np.zeros(self.nlocal_parameters) # Mean values of local parameters
        self.B = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Cholesky decomposition of covariance matrices
        self.Bi = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Inverses of the Cholesky decomposition of covariance matrices
        self.prior_product = np.zeros(1)
        self.ppratio = 0
        
        # RJMCMC
        self.likelihood = 0.
        
        # Global function values
        self.gf_values = np.zeros(self.npoints)
        self.global_function = global_function
        self.lf_all_domain = np.zeros(self.npoints)
        
        # For internal computations
        self._u = np.zeros(self.nlocal_parameters) # Random vector
        self._mu = np.zeros(self.nlocal_parameters) # Mean of local parameters
        self._sigma = np.ones(self.nlocal_parameters) # Standard deviation of local parameters
        self._S = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Covariance matrices of each local model
        self._B = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Cholesky decomposition of covariance matrix
        self._Bi = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Inverse of the Cholesky decomposition of covariance matrix
        self._N = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Normal matrix
        self._Q = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Cofactor matrix (N^T*N)^-1
        self._L = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Lower triangular matrix
        self._U = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Upper triangular matrix
        self._P = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Permutation matrix
        self._aux1 = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Auxiliar matrix
        self._aux2 = np.zeros((self.nlocal_parameters, self.nlocal_parameters)) # Auxiliar matrix
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
            self.X = self.design_matrices(dataset[:,0], self.theta_g)
        else:
            self.X = self.design_matrices.copy()
        
        # Compute model for each segment
        self._u = self.nprs.normal(size=self.nlocal_parameters)
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
        self._u = np.dot(self.Bi, self.theta - self.mu)
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
            self.X = self.design_matrices(dataset[:,0], self.theta_g)
        
        prob = 1./self.ppratio
        self._u = np.dot(self.Bi, self.theta - self.mu)
        
        curve_prob = self._update_partition(dataset)
        if curve_prob <= 0:
            return 0, 0.
            
        prob *= curve_prob
        
        return 1, prob
    
    def misfit(self, datasets):
        
        dataset = datasets[0]
        
        if self.nglobal_parameters > 0:
            likelihood = self.likelihood_function(self.theta_g, self.theta, dataset[:,0], dataset[:,1], dataset[:,2])
        else:
            likelihood = self.likelihood_function(self.theta, dataset[:,0], dataset[:,1], dataset[:,2])
            
        self.likelihood = likelihood
        
    def sample(self, x, di):
        
        if self.nglobal_parameters > 0:
            return self.sampling_function(x, self.theta_g, self.theta)
        else:
            return self.sampling_function(x, self.theta)
    
    def _update_partition(self, dataset):
        # Compute OLS and determine model probabilities
        prob = utils.update_partition_single_model(dataset,
                                                   self.gf_values,
                                                   self.lf_all_domain,
                                                   0,
                                                   self.npoints,
                                                   self.nlocal_parameters,
                                                   self._mu,
                                                   self._sigma,
                                                   self._S,
                                                   self.X,
                                                   self._N,
                                                   self._Q,
                                                   self._L,
                                                   self._U,
                                                   self._P,
                                                   self._aux1,
                                                   self._aux2,
                                                   self._auto_z,
                                                   self.prior_product,
                                                   self.B,
                                                   self.Bi,
                                                   self._u,
                                                   self.mu,
                                                   self.theta)
        
        self.ppratio = prob
        
        return prob
