# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 15:47:22 2021

@author: sebac
"""

import numpy as np

class mp_so_iu:
    
    def __init__(self, npoints, ndatasets, xmin, xmax, min_partitions, max_partitions,
                 kmax, samples, pd, options, weights, nglobal_parameters,
                 global_parameters, nlocal_parameters, local_parameters,
                 computeDesignMatrices, loglikelihood_function, global_function,
                 design_matrices, local_function, sampling_function, fixed_init, nprs):
        
        # Datasets parameters
        self.npoints = npoints.copy()
        self.ndatasets = ndatasets
        self.xmin = xmin
        self.xmax = xmax
        
        # Models parameters
        self.nlocal_parameters = nlocal_parameters
        self.local_parameters = local_parameters.copy()
        
        self.nglobal_parameters = nglobal_parameters
        if self.nglobal_parameters > 0:
            self.global_parameters = global_parameters.copy()
        
        # Partition parameters
        self.npartitions = 0
        self.min_partitions = min_partitions
        self.max_partitions = max_partitions
        self.partition_map = [0]
        
        # For the sake of simplicity, we define boundaries
        self.nboundaries = 0
        self.boundaries = []
        self.min_boundaries = min_partitions + 1
        self.max_boundaries = max_partitions + 1
        
        # Boundary perturbation
        self.pd = pd
        
        # Model values
        self.theta_g = np.empty((self.ndatasets, self.nglobal_parameters))
        self.theta = np.zeros((self.ndatasets, self.max_partitions, self.nlocal_parameters)) # Sampled values of local parameters
        
        # RJMCMC
        self.likelihood = 0.
        
        self.nprs = nprs
        
        # Internal calculations
        self.delta_prod = (self.local_parameters[:,:,1] - self.local_parameters[:,:,0]).prod()
        
        # Functions
        self.sampling_function = sampling_function
        self.likelihood_function = loglikelihood_function
    
    def initialize(self, datasets):
        # By default the first model uses the minimum partitions allowed
        self.npartitions = self.min_partitions
        self.nboundaries = self.min_boundaries
        self.boundaries = [self.xmin, self.xmax]
        
        # Create internal boundaries if necessary
        bi = 2
        while bi < self.nboundaries:
            # Sample random boundary location
            x = self.nprs.random()*(self.xmax - self.xmin) + self.xmin
            
            # Indentify position with other boundaries
            i = 0
            while x > self.boundaries[i]:
                i += 1
                
            # Check if there is data available in the new proposed partitions
            if self._check_data(i, x, datasets):
                self.boundaries.insert(i, x)
                self.partition_map.insert(i, bi-1)
                bi += 1
        
        # Sample global parameters if there is any
        if self.nglobal_parameters > 0:
            self.theta_g = self.nprs.random((self.ndatasets, self.nglobal_parameters))*\
                (self.global_parameters[:,:,1] - self.global_parameters[:,:,0]) + self.global_parameters[:,:,0]
        
        # Sample local parameters
        self.theta = self.nprs.random((self.ndatasets, self.max_partitions, self.nlocal_parameters))*\
            (self.local_parameters[:,np.newaxis,:,1] - self.local_parameters[:,np.newaxis,:,0]) + self.local_parameters[:,np.newaxis,:,0]

    def initialize_fixed(self, datasets, fixed_init):
        
        if 'local_parameters' in fixed_init:
            #TODO: check size
            self.theta = np.array(fixed_init['local_parameters'])
        else:
            self.theta = self.nprs.random((self.ndatasets, self.max_partitions, self.nlocal_parameters))*\
                (self.local_parameters[:,np.newaxis,:,1] - self.local_parameters[:,np.newaxis,:,0]) + self.local_parameters[:,np.newaxis,:,0]
            
        return 1
            
    def clone(self, other):
        
        # Partition parameters
        self.npartitions = other.npartitions
        self.min_partitions = other.min_partitions
        self.max_partitions = other.max_partitions
        self.partition_map = other.partition_map.copy()
        
        # For the sake of simplicity, we define boundaries
        self.nboundaries = other.nboundaries
        self.boundaries = other.boundaries.copy()
        self.min_boundaries = other.min_boundaries
        self.max_boundaries = other.max_boundaries
        
        # Boundary perturbation
        self.pd = other.pd
        
        # Model values
        self.nglobal_parameters = other.nglobal_parameters
        if self.nglobal_parameters > 0:
            self.theta_g = other.theta_g.copy()
        
        self.nlocal_parameters = other.nlocal_parameters
        self.theta = other.theta.copy() #[other.theta[di].copy() for di in range(other.ndatasets)]
        
        # RJMCMC
        self.likelihood = other.likelihood
        
    
    def perturb(self, process, datasets):
        
        if process == 0:
            status, process_prob = self._propose_birth(datasets)
        elif process == 1:
            status, process_prob = self._propose_death(datasets)
        elif process == 2:
            status, process_prob = self._propose_move(datasets)
        elif process == 3:
            status, process_prob = self._propose_merge(datasets)
        elif process == 4:
            status, process_prob = self._propose_local_parameters(datasets)
        else:
            status, process_prob = self._propose_global_parameters(datasets)
        
        return status, process_prob
    
    def _propose_birth(self, datasets):
        
        if self.npartitions == self.max_partitions:
            return 0, 0.
        
        new_x = self.nprs.random()*(self.xmax - self.xmin) + self.xmin
        new_iy = 1
        while new_iy in self.partition_map:
            new_iy += 1
        
        i = 0
        while new_x > self.boundaries[i]:
            i += 1
        # prev_i = i-1
        
        if not self._check_data(i, new_x, datasets):
            return 0, 0.
        
        prob = 1.
            
        self.boundaries.insert(i, new_x)
        self.partition_map.insert(i, new_iy)
        self.npartitions += 1
        self.nboundaries += 1
        
        self.theta[:,new_iy] = self.nprs.random((self.ndatasets, self.nlocal_parameters))*\
            (self.local_parameters[:,:,1] - self.local_parameters[:,:,0]) + self.local_parameters[:,:,0]
            
        prob = 1.
        
        return 1, prob
    
    def _propose_death(self, datasets):
        
        if self.npartitions <= self.min_partitions:
            return 0, 0.
        
        del_iy = self.nprs.choice(self.npartitions-1) + 1
                
        self.partition_map.pop(del_iy)
        self.boundaries.pop(del_iy)
        self.npartitions -= 1
        self.nboundaries -= 1
        
        prob = 1.
            
        return 1, prob
    
    def _propose_move(self, datasets):
        
        if self.npartitions <= self.min_partitions:
            return 0, 0.
        
        iy = self.nprs.choice(self.npartitions-1) + 1
        
        old_x = self.boundaries[iy]
        new_x = old_x + self.nprs.normal(0., self.pd)
        
        xi = self.boundaries[iy-1]
        xf = self.boundaries[iy+1]
        
        if new_x <= xi or new_x >= xf:
            return 0, 0.
        
        for dataset in datasets:
            n = len(np.where((dataset[:,0] >= xi) & (dataset[:,0] < new_x))[0])
            if n <= self.nlocal_parameters:
                return 0, 0.
            n = len(np.where((dataset[:,0] >= new_x) & (dataset[:,0] <= xf))[0])
            if n <= self.nlocal_parameters:
                return 0, 0.
        
        self.boundaries[iy] = new_x
        
        prob = 1.
            
        return 1, prob
    
    def _propose_merge(self, datasets):
        
        if self.npartitions <= self.min_partitions or self.nboundaries <= 3:
            return 0, 0.
        
        jy = self.nprs.choice(self.npartitions-1) + 1
        
        iyop = jy - 1
        iynp = jy + 1
        
        iy = iyop
        ky = iynp
        if iyop == 0:
            case = 1
        elif iynp == self.npartitions:
            case = 2
        else:
            case = 3
        
        xa = self.boundaries[jy]
        if case == 1:
            ik = 1
            xb = self.boundaries[ky]
            # prev_i = iy
        elif case == 2:
            ik = 0
            xb = self.boundaries[iy]
            # prev_i = iy - 1
        else:
            ik = self.nprs.integers(0, 2)
            if ik == 0:
                xb = self.boundaries[iy]
                # prev_i = iy - 1
            else:
                xb = self.boundaries[ky]
                # prev_i = iy
        
        alpha = self.nprs.random()
        new_x = xa*alpha + xb*(1.-alpha)
        
        if ik == 0:
            self.boundaries[iy] = new_x
        else:
            self.boundaries[ky] = new_x
            
        self.partition_map.pop(jy)
        self.boundaries.pop(jy)
        self.npartitions -= 1
        self.nboundaries -= 1
        
        prob = 1.
            
        return 1, prob
    
    def _propose_local_parameters(self, datasets):
        
        if self.ndatasets == 1:
            di = 0
        else:
            di = self.nprs.choice(self.ndatasets)
        
        if self.nlocal_parameters == 1:
            li = 0
        else:
            li = self.nprs.choice(self.nlocal_parameters)
            
        pi = self.nprs.choice(self.npartitions)
        
        self.theta[di,pi,li] = self.nprs.random()*(self.local_parameters[di,li,1] - self.local_parameters[di,li,0])\
            + self.local_parameters[di,li,0] 
        
        prob = 1.
        
        return 1, prob
    
    def _propose_global_parameters(self, datasets):
        
        if self.ndatasets == 1:
            di = 0
        else:
            di = self.nprs.choice(self.ndatasets)
        
        if self.nglobal_parameters == 1:
            gi = 0
        else:
            gi = self.nprs.choice(self.nglobal_parameters)
        
        self.theta_g[di][gi] = self.nprs.random()*(self.global_parameters[di,gi,1] - self.global_parameters[di,gi,0])\
            + self.global_parameters[di,gi,0]
        
        prob = 1.
        
        return 1, prob
    
    def _check_data(self, i, x, datasets):
        xi = self.boundaries[i-1]
        xf = self.boundaries[i]
        
        for dataset in datasets:
            n = len(np.where((dataset[:,0] >= xi) & (dataset[:,0] < x))[0])
            if n <= self.nlocal_parameters:
                return False
            n = len(np.where((dataset[:,0] >= x) & (dataset[:,0] <= xf))[0])
            if n <= self.nlocal_parameters:
                return False
        
        return True

    def misfit(self, datasets):
        
        likelihood = 0.
        for di, dataset in enumerate(datasets):
            if self.nglobal_parameters > 0:
                likelihood += self.likelihood_function(self.theta_g[di],
                                                       self.boundaries,
                                                       self.theta[di,self.partition_map],
                                                       dataset[:,0],
                                                       dataset[:,1],
                                                       dataset[:,2])
            else:
                likelihood += self.likelihood_function(self.boundaries,
                                                       self.theta[di,self.partition_map],
                                                       dataset[:,0],
                                                       dataset[:,1],
                                                       dataset[:,2])
            
        self.likelihood = likelihood
    
    def sample(self, x, di):
        if self.nglobal_parameters > 0:
            return self.sampling_function(x, self.theta_g[di], self.boundaries, self.theta[di,self.partition_map])
        else:
            return self.sampling_function(x, self.boundaries, self.theta[di,self.partition_map])
