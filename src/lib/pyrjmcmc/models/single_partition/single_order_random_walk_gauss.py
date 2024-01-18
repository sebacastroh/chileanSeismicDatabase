# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 21:53:34 2021

@author: sebac
"""

import numpy as np

class sp_so_rwg:
    
    def __init__(self, npoints, ndatasets, xmin, xmax, min_partitions, max_partitions,
                 kmax, samples, pd, options, weights, nglobal_parameters,
                 global_parameters, nlocal_parameters, local_parameters,
                 computeDesignMatrices, loglikelihood_function, global_function,
                 design_matrices, local_function, sampling_function, fixed_init, nprs):
        

        # Datasets parameters
        self.npoints = npoints[0]
        
        # Models parameters
        self.nlocal_parameters = nlocal_parameters
        self.local_parameters = local_parameters.copy()
        
        self.theta = np.empty(self.nlocal_parameters)
        
        # RJMCMC
        self.likelihood = 0.
        
        self.nprs = nprs
        
        # Functions
        self.sampling_function = sampling_function
        self.likelihood_function = loglikelihood_function
    
    def initialize(self, datasets):
        # Sample parameters
        self.theta = self.nprs.random(self.nlocal_parameters)*\
            (self.local_parameters[:,1] - self.local_parameters[:,0]) + self.local_parameters[:,0]

    def initialize_fixed(self, datasets, fixed_init):
        
        if 'local_parameters' in fixed_init:
            #TODO: check size
            self.theta = np.array(fixed_init['local_parameters'])
        else:
            self.theta = self.nprs.random((self.ndatasets, self.nlocal_parameters))*\
                (self.local_parameters[:,np.newaxis,:,1] - self.local_parameters[:,np.newaxis,:,0]) + self.local_parameters[:,np.newaxis,:,0]
            
        return 1
            
    def clone(self, other):
        
        # Model values
        self.theta = other.theta.copy()
        
        # RJMCMC
        self.likelihood = other.likelihood
        
    
    def perturb(self, process, datasets):
        
        dataset = datasets[0]
        status, process_prob = self._propose_parameters(dataset)
        
        return status, process_prob
    
    def _propose_parameters(self, dataset):
        
        if self.nlocal_parameters == 1:
            li = 0
        else:
            li = self.nprs.integers(0, self.nlocal_parameters)
        
        dv = self.nprs.normal(0., self.local_parameters[li,2])
        
        self.theta[li] += dv
        
        if self.theta[li] < self.local_parameters[li,0] or\
           self.theta[li] > self.local_parameters[li,1]:
               
           return 0, 0. 

        prob = 1.
        
        return 1, prob

    def misfit(self, datasets):
        
        dataset = datasets[0]
        likelihood = self.likelihood_function(self.theta, dataset[:,0], dataset[:,1], dataset[:,2])
            
        self.likelihood = likelihood
    
    def sample(self, x, di):
        return self.sampling_function(x, self.theta)
