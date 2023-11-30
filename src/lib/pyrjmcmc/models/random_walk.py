# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 16:10:27 2021

@author: sebac
"""

import numpy as np

class single_order_custom_likelihood_random_walk:
    
    def __init__(self, npoints, ndatasets,
                 nglobal_parameters, nlocal_parameters,
                 min_partitions, max_partitions,
                 xmin, xmax, pd, nprs, rs):

        # Datasets parameters
        self.npoints = npoints.copy()
        self.ndatasets = ndatasets
        self.xmin = xmin
        self.xmax = xmax

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
        
        # Models parameters
        self.nglobal_parameters = nglobal_parameters
        self.nlocal_parameters = nlocal_parameters
        self.global_parameters = np.empty((self.ndatasets, self.nglobal_parameters))
        self.theta = np.empty((self.ndatasets, self.max_partitions, self.nlocal_parameters))
        
        # Boundary perturbation
        self.pd = pd
        
        # RJMCMC
        self.likelihood = 0.
        
        self.nprs = nprs
        self.rs = rs
        
        self.k = np.zeros((self.ndatasets, self.max_partitions), dtype=int)
    
    def initialize(self, datasets, global_parameters, local_parameters):
        # By default the first model uses the minimum partitions allowed
        self.npartitions = self.min_partitions
        self.nboundaries = self.min_boundaries
        self.boundaries = [self.xmin, self.xmax]
        
        # Create internal boundaries if necessary
        bi = 2
        while bi < self.nboundaries:
            # Sample random boundary location
            x = self.rs.random()*(self.xmax - self.xmin) + self.xmin
            
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
            self.global_parameters = self.nprs.rand(self.ndatasets, self.nglobal_parameters)*\
                (global_parameters[:,:,1] - global_parameters[:,:,0]) + global_parameters[:,:,0]
        
        # Sample local parameters
        self.theta[:,:self.npartitions,:] = self.nprs.rand(self.ndatasets, self.npartitions, self.nlocal_parameters)*\
            (local_parameters[:,np.newaxis,:,1] - local_parameters[:,np.newaxis,:,0]) + local_parameters[:,np.newaxis,:,0]

    def initialize_fixed(self, datasets, fixed_init, global_parameters, local_parameters):
        
        if 'boundaries' in fixed_init:
            boundaries = fixed_init['boundaries']
        else:
            boundaries = [self.xmin, self.xmax]
        
        i = 0
        j = len(boundaries)
        if boundaries[0] == self.xmin:
            i = 1
        if boundaries[-1] == self.xmax:
            j -= 1
        
        self.boundaries = [self.xmin]
        self.partition_map = [0]
        for bi, b in enumerate(boundaries[i:j]):
            self.boundaries.append(b)
            self.partition_map.append(bi+1)
        self.boundaries.append(self.xmax)
        
        self.npartitions = len(self.boundaries) - 1
        self.nboundaries = len(self.boundaries)
        
        if self.nglobal_parameters > 0:
            if 'global_parameters' in fixed_init:
                self.global_parameters = np.array(fixed_init['global_parameters'])
            else:
                self.global_parameters = self.nprs.rand(self.ndatasets, self.nglobal_parameters)*\
                    (global_parameters[:,:,1] - global_parameters[:,:,0]) + global_parameters[:,:,0]
        

        if 'local_parameters' in fixed_init:
            #TODO: check size
            self.theta = np.array(fixed_init['local_parameters'])
        else:
            self.theta = self.nprs.rand(self.ndatasets, self.nlocal_parameters)*\
                (local_parameters[:,np.newaxis,:,1] - local_parameters[:,np.newaxis,:,0]) + local_parameters[:,np.newaxis,:,0]
            
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
        self.global_parameters = other.global_parameters.copy()
        self.theta = other.theta.copy()
        
        # RJMCMC
        self.likelihood = other.likelihood
        
    
    def perturb(self, process, datasets, global_parameters, local_parameters):
        
        if process == 0:
            status, process_prob = self.propose_birth(datasets, local_parameters)
        elif process == 1:
            status, process_prob = self.propose_death(local_parameters)
        elif process == 2:
            status, process_prob = self.propose_move(datasets, local_parameters)
        elif process == 3:
            status, process_prob = self.propose_merge(datasets, local_parameters)
        elif process == 4:
            status, process_prob = self.propose_local_parameters(datasets, local_parameters)
        else:
            status, process_prob = self.propose_global_parameters(datasets, global_parameters)
        
        return status, process_prob
    
    def propose_birth(self, datasets, local_parameters):
        
        if self.npartitions == self.max_partitions:
            return 0, 0.
        
        new_x = self.rs.random()*(self.xmax - self.xmin) + self.xmin
        new_iy = 1
        while new_iy in self.partition_map:
            new_iy += 1
        
        i = 0
        while new_x > self.boundaries[i]:
            i += 1
        prev_i = self.partition_map[i-1]
        
        if not self._check_data(i, new_x, datasets):
            return 0, 0.
        
        self.boundaries.insert(i, new_x)
        self.partition_map.insert(i, new_iy)
        self.npartitions += 1
        self.nboundaries += 1
        
        dv = self.nprs.normal(0., local_parameters[:,:,2])
        self.theta[:,new_iy] = self.theta[:,prev_i] + dv

        if np.any(self.theta[:,new_iy] < local_parameters[:,:,0]):
            return 0, 0.

        if np.any(self.theta[:,new_iy] > local_parameters[:,:,1]):
            return 0, 0.

        prob = (np.exp(-0.5*(dv/local_parameters[:,:,2])**2)/(local_parameters[:,:,2]*np.sqrt(2.*np.pi))).prod()
        prob = 1./prob
        if prob == 0.:
            return 0, 0.

        return 1, prob
    
    def propose_death(self, local_parameters):
        
        if self.npartitions <= self.min_partitions:
            return 0, 0.
        
        del_iy = self.nprs.choice(self.npartitions-1) + 1
        
        prev_i = del_iy-1

        deleted_parameters = self.theta[:,self.partition_map[del_iy],:]
                
        self.partition_map.pop(del_iy)
        self.boundaries.pop(del_iy)
        self.npartitions -= 1
        self.nboundaries -= 1
                
        prob = 1./(local_parameters[:,:,1] - local_parameters[:,:,0]).prod()

        dv = self.theta[:,self.partition_map[prev_i]] - deleted_parameters
        prob = (np.exp(-0.5*(dv/local_parameters[:,:,2])**2)/(local_parameters[:,:,2]*np.sqrt(2.*np.pi))).prod()
        if prob == 0.:
            return 0, 0.

        return 1, prob
    
    def propose_move(self, datasets, local_parameters):
        
        if self.npartitions <= self.min_partitions:
            return 0, 0.
        
        iy = self.nprs.choice(self.npartitions-1) + 1
        
        old_x = self.boundaries[iy]
        new_x = old_x + self.rs.gauss(0., self.pd)
        
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
    
    def propose_merge(self, datasets, local_parameters):
        
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
            prev_i = iy
        elif case == 2:
            ik = 0
            xb = self.boundaries[iy]
            prev_i = iy - 1
        else:
            ik = self.rs.randint(0, 1)
            if ik == 0:
                xb = self.boundaries[iy]
                prev_i = iy - 1
            else:
                xb = self.boundaries[ky]
                prev_i = iy
        
        alpha = self.rs.random()
        new_x = xa*alpha + xb*(1.-alpha)
        
        if ik == 0:
            self.boundaries[iy] = new_x
        else:
            self.boundaries[ky] = new_x
        
        deleted_parameters = self.theta[:,self.partition_map[jy],:]    
        
        self.partition_map.pop(jy)
        self.boundaries.pop(jy)
        self.npartitions -= 1
        self.nboundaries -= 1
        
        dv = self.theta[:,self.partition_map[prev_i]] - deleted_parameters
        prob = (np.exp(-0.5*(dv/local_parameters[:,:,2])**2)/(local_parameters[:,:,2]*np.sqrt(2.*np.pi))).prod()
        
        if prob == 0.:
            return 0, 0.
        
        return 1, prob
    
    def propose_local_parameters(self, datasets, local_parameters):
        
        if self.ndatasets == 1:
            di = 0
        else:
            di = self.nprs.choice(self.ndatasets)
            
        iy = self.nprs.choice(self.npartitions)

        if self.nlocal_parameters == 0:
            li = 0
        else:
            li = self.nprs.choice(self.nlocal_parameters)

        dv = self.rs.gauss(0., local_parameters[di,li,2])
        self.theta[di, self.partition_map[iy], li] += dv
        
        if self.theta[di, self.partition_map[iy], li] < local_parameters[di,li,0]:
            return 0, 0.

        if self.theta[di, self.partition_map[iy], li] > local_parameters[di,li,1]:
            return 0, 0.

        prob = 1.
        
        return 1, prob
    
    def propose_global_parameters(self, datasets, global_parameters):
        
        if self.ndatasets == 1:
            di = 0
        else:
            di = self.nprs.choice(self.ndatasets)
            
        if self.nglobal_parameters == 1:
            gi = 0
        else:
            gi = self.nprs.choice(self.nglobal_parameters)
            
        self.global_parameters[di,gi] += self.rs.gauss(0., global_parameters[di,gi,2])
        
        if self.global_parameters[di,gi] < global_parameters[di,gi,0] or \
            self.global_parameters[di,gi] > global_parameters[di,gi,1]:
            return 0, 0.
        
        prob = 1.
        
        return 1, prob

    def misfit(self, datasets: list, likelihood_function):
        
        likelihood = 0.
        for di, dataset in enumerate(datasets):
            likelihood += likelihood_function(self.global_parameters[di],
                                              self.boundaries,
                                              self.theta[di][self.partition_map],
                                              dataset[:,0],
                                              dataset[:,1],
                                              dataset[:,2])
            
        self.likelihood = likelihood
    
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
    