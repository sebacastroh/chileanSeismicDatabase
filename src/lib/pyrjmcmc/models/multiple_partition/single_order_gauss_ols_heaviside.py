# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 11:57:45 2021

@author: sebac
"""

import utils
import numpy as np

class mp_so_goh:
    
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
        self.nglobal_parameters = nglobal_parameters
        self.nlocal_parameters = nlocal_parameters
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
        
        # Global function values
        self.gf_values = [np.zeros(npoints[di]) for di in range(ndatasets)]
        
        # Design matrices
        self.computeDesignMatrices = computeDesignMatrices
        self.X = [np.zeros((npoints[di], self.nlocal_parameters)) for di in range(ndatasets)]
        self.design_matrices = design_matrices
        
        # Local functions values
        self.lf_all_domain = [np.zeros(npoints[di]) for di in range(ndatasets)]
        self.lf_values = [np.zeros((npoints[di], self.max_partitions)) for di in range(ndatasets)]
        
        # Current local models
        self.theta = [np.zeros((self.max_partitions, self.nlocal_parameters)) for di in range(ndatasets)] # Sampled values of local parameters
        self.mu = [np.zeros((self.max_partitions, self.nlocal_parameters)) for di in range(ndatasets)] # Mean values of local parameters
        self.B = [np.zeros((self.max_partitions, self.nlocal_parameters, self.nlocal_parameters)) for di in range(ndatasets)] # Cholesky decomposition of covariance matrices
        self.Bi = [np.zeros((self.max_partitions, self.nlocal_parameters, self.nlocal_parameters)) for di in range(ndatasets)] # Inverses of the Cholesky decomposition of covariance matrices
        self.prior_product = [np.zeros(self.max_partitions) for di in range(ndatasets)]
        self.ppratio = [np.zeros(self.max_partitions) for di in range(ndatasets)]
        
        # RJMCMC
        self.likelihood = 0.
        
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
        
        self.nprs = nprs
        
        # Functions
        self.sampling_function = sampling_function
        self.likelihood_function = loglikelihood_function
        self.global_function = global_function
        self.local_function = local_function
    
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
            
            #Evaluate global function
            if callable(self.global_function):
                for di, dataset in enumerate(datasets):
                    self.gf_values[di] = self.global_function(dataset[:,0], self.theta_g[di])
        
        # Compute design matrices if corresponds
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                self.X[di] = self.design_matrices(dataset[:,0], self.theta_g[di], self.boundaries)
        else:
            for di in range(self.ndatasets):
                self.X[di] = self.design_matrices[di].copy()
        
        # Compute model for each segment
        for i, pi in enumerate(self.partition_map):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            
            for di, dataset in enumerate(datasets):
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                self._u = self.nprs.normal(size=self.nlocal_parameters)
                
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                self._update_partition(datasets, di, pi, xi, xj)
                
                if self.nglobal_parameters == 0:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], xl, xr, self.theta[di][pi])
                else:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], self.theta_g[di], xl, xr, self.theta[di][pi])
                
        
    def initialize_fixed(self, datasets, fixed_init):
        
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
                self.theta_g = np.array(fixed_init['global_parameters'])
            else:
                self.theta_g = self.nprs.random((self.ndatasets, self.nglobal_parameters))*\
                    (self.global_parameters[:,:,1] - self.global_parameters[:,:,0]) + self.global_parameters[:,:,0]
                
            for di, dataset in enumerate(datasets):
                self.gf_values[di] = self.global_function(dataset[:,0], self.theta_g[di])
            
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                self.X[di] = self.design_matrices(dataset[:,0], self.theta_g[di], self.boundaries)
        else:
            for di in range(self.ndatasets):
                self.X[di] = self.design_matrices[di].copy()
                
        # Compute model for each segment
        for i, pi in enumerate(self.partition_map):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            
            for di, dataset in enumerate(datasets):
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                self._u = self.nprs.normal(size=self.nlocal_parameters)
                
                if 'local_parameters' in fixed_init and 'boundaries' in fixed_init:
                    self.theta[di][pi] = fixed_init['local_parameters'][di][pi]
                    use_theta = 1
                    prob = self._update_partition_fixed(datasets, di, pi, xi, xj, use_theta)
                else:
                    prob = self._update_partition(datasets, di, pi, xi, xj)
                
                if prob <= 0.:
                    return 0
            
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
        
        # Global function values
        if self.nglobal_parameters > 0:
            self.gf_values = other.gf_values.copy()
        
        # Design matrices
        self.computeDesignMatrices = other.computeDesignMatrices
        self.X = [other.X[di].copy() for di in range(other.ndatasets)]
        
        # Local functions values
        self.lf_all_domain = [other.lf_all_domain[di].copy() for di in range(other.ndatasets)]
        self.lf_values = [other.lf_values[di].copy() for di in range(other.ndatasets)]
        
        # Current local models
        self.theta = [other.theta[di].copy() for di in range(other.ndatasets)]
        self.mu = [other.mu[di].copy() for di in range(other.ndatasets)]
        self.B = [other.B[di].copy() for di in range(other.ndatasets)]
        self.Bi = [other.Bi[di].copy() for di in range(other.ndatasets)]
        self.prior_product = [other.prior_product[di].copy() for di in range(other.ndatasets)]
        self.ppratio = [other.ppratio[di].copy() for di in range(other.ndatasets)]
        
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
        prev_i = i-1
        
        if not self._check_data(i, new_x, datasets):
            return 0, 0.
        
        prob = 1.
        
        # Save prior product
        for di in range(self.ndatasets):
            prob /= self.ppratio[di][self.partition_map[prev_i:]].prod()
        
        self.boundaries.insert(i, new_x)
        self.partition_map.insert(i, new_iy)
        self.npartitions += 1
        self.nboundaries += 1
        
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                self.X[di] = self.design_matrices(dataset[:,0], self.theta_g[di], self.boundaries)
        
        for i, pi in enumerate(self.partition_map[prev_i:], prev_i):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            for di, dataset in enumerate(datasets):
                
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                
                if pi == new_iy:
                    ppi = self.partition_map[prev_i]
                    self._u = np.dot(self.Bi[di][ppi], self.theta[di][ppi] - self.mu[di][ppi])
                else:
                    self._u = np.dot(self.Bi[di][pi], self.theta[di][pi] - self.mu[di][pi])
                
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                curve_prob = self._update_partition(datasets, di, pi, xi, xj)
                
                if curve_prob <= 0:
                    return 0, 0.
                
                prob *= curve_prob
                
                if self.nglobal_parameters == 0:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], xl, xr, self.theta[di][pi])
                else:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], self.theta_g[di], xl, xr, self.theta[di][pi])
        
        return 1, prob
                
    
    def _propose_death(self, datasets):
        
        if self.npartitions <= self.min_partitions:
            return 0, 0.
        
        del_iy = self.nprs.choice(self.npartitions-1) + 1
        prev_i = del_iy-1
        
        prob = 1.
        
        for di in range(self.ndatasets):
            prob /= self.ppratio[di][self.partition_map[prev_i:]].prod()
                
        self.partition_map.pop(del_iy)
        self.boundaries.pop(del_iy)
        self.npartitions -= 1
        self.nboundaries -= 1
        
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                self.X[di] = self.design_matrices(dataset[:,0], self.theta_g[di], self.boundaries)
        
        for i, pi in enumerate(self.partition_map[prev_i:], prev_i):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            for di, dataset in enumerate(datasets):
                
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                
                self._u = np.dot(self.Bi[di][pi], self.theta[di][pi] - self.mu[di][pi])
                    
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                curve_prob = self._update_partition(datasets, di, pi, xi, xj)
                
                if curve_prob <= 0:
                    return 0, 0.
                
                prob *= curve_prob
                
                if self.nglobal_parameters == 0:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], xl, xr, self.theta[di][pi])
                else:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], self.theta_g[di], xl, xr, self.theta[di][pi])
            
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
        
        prev_i = iy - 1
        self.boundaries[iy] = new_x
        
        prob = 1.
        
        for di in range(self.ndatasets):
            prob /= self.ppratio[di][self.partition_map[prev_i:]].prod()
            
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                self.X[di] = self.design_matrices(dataset[:,0], self.theta_g[di], self.boundaries)
        
        for i, pi in enumerate(self.partition_map[prev_i:], prev_i):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            for di, dataset in enumerate(datasets):
                
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                
                self._u = np.dot(self.Bi[di][pi], self.theta[di][pi] - self.mu[di][pi])
                    
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                curve_prob = self._update_partition(datasets, di, pi, xi, xj)
                
                if curve_prob <= 0:
                    return 0, 0.
                
                prob *= curve_prob
                
                if self.nglobal_parameters == 0:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], xl, xr, self.theta[di][pi])
                else:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], self.theta_g[di], xl, xr, self.theta[di][pi])
                
            
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
            prev_i = iy
        elif case == 2:
            ik = 0
            xb = self.boundaries[iy]
            prev_i = iy - 1
        else:
            ik = self.nprs.integers(0, 2)
            if ik == 0:
                xb = self.boundaries[iy]
                prev_i = iy - 1
            else:
                xb = self.boundaries[ky]
                prev_i = iy
        
        alpha = self.nprs.random()
        new_x = xa*alpha + xb*(1.-alpha)
        
        prob = 1.
        for di in range(self.ndatasets):
            prob /= self.ppratio[di][self.partition_map[prev_i:]].prod()
        
        if ik == 0:
            self.boundaries[iy] = new_x
        else:
            self.boundaries[ky] = new_x
            
        self.partition_map.pop(jy)
        self.boundaries.pop(jy)
        self.npartitions -= 1
        self.nboundaries -= 1
        
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                self.X[di] = self.design_matrices(dataset[:,0], self.theta_g[di], self.boundaries)
        
        for i, pi in enumerate(self.partition_map[prev_i:], prev_i):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            for di, dataset in enumerate(datasets):
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                
                self._u = np.dot(self.Bi[di][pi], self.theta[di][pi] - self.mu[di][pi])
                    
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                curve_prob = self._update_partition(datasets, di, pi, xi, xj)
                
                if curve_prob <= 0:
                    return 0, 0.
                
                prob *= curve_prob
                
                if self.nglobal_parameters == 0:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], xl, xr, self.theta[di][pi])
                else:
                    self.lf_values[di][:,pi] = self.local_function(dataset[:,0], self.theta_g[di], xl, xr, self.theta[di][pi])
        
        
        return 1, prob
    
    def _propose_local_parameters(self, datasets):
        
        if self.ndatasets == 1:
            di = 0
        else:
            di = self.nprs.choice(self.ndatasets)
            
        iy = self.nprs.choice(self.npartitions)
        prob = 1./self.ppratio[di][self.partition_map[iy:]].prod()
        
        for i, pi in enumerate(self.partition_map[iy:], iy):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            dataset = datasets[di]
            dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
            xi, xj = dataset_indices[[0, -1]]
            
            self._u = np.dot(self.Bi[di][pi], self.theta[di][pi] - self.mu[di][pi])
                
            self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
            curve_prob = self._update_partition(datasets, di, pi, xi, xj)
            
            if curve_prob <= 0:
                return 0, 0.
            
            prob *= curve_prob
            
            if self.nglobal_parameters == 0:
                self.lf_values[di][:,pi] = self.local_function(dataset[:,0], xl, xr, self.theta[di][pi])
            else:
                self.lf_values[di][:,pi] = self.local_function(dataset[:,0], self.theta_g[di], xl, xr, self.theta[di][pi])
        
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
            
        self.theta_g[di,gi] += self.nprs.normal(0., self.global_parameters[di,gi,2])
        
        if self.theta_g[di,gi] < self.global_parameters[di,gi,0] or \
            self.theta_g[di,gi] > self.global_parameters[di,gi,1]:
            return 0, 0.
        
        if callable(self.global_function):
            self.gf_values[di] = self.global_function(datasets[di][:,0], self.theta_g[di])
        
        if self.computeDesignMatrices:
            self.X[di] = self.design_matrices(datasets[di][:,0], self.theta_g[di], self.boundaries)
        
        prob = 1./self.ppratio[di][self.partition_map].prod()
        
        for i, pi in enumerate(self.partition_map):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            dataset = datasets[di]
            dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
            xi, xj = dataset_indices[[0, -1]]
            
            self._u = np.dot(self.Bi[di][pi], self.theta[di][pi] - self.mu[di][pi])
                
            self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
            curve_prob = self._update_partition(datasets, di, pi, xi, xj)
            
            if curve_prob <= 0:
                return 0, 0.
            
            prob *= curve_prob
            
            if self.nglobal_parameters == 0:
                self.lf_values[di][:,pi] = self.local_function(dataset[:,0], xl, xr, self.theta[di][pi])
            else:
                self.lf_values[di][:,pi] = self.local_function(dataset[:,0], self.theta_g[di], xl, xr, self.theta[di][pi])
            
        return 1, prob
    
    def misfit(self, datasets):
        
        likelihood = 0.
        for di, dataset in enumerate(datasets):
            if self.nglobal_parameters > 0:
                likelihood += self.likelihood_function(self.theta_g[di],
                                                       self.boundaries,
                                                       self.theta[di][self.partition_map],
                                                       dataset[:,0],
                                                       dataset[:,1],
                                                       dataset[:,2])
            else:
                likelihood += self.likelihood_function(self.boundaries,
                                                       self.theta[di][self.partition_map],
                                                       dataset[:,0],
                                                       dataset[:,1],
                                                       dataset[:,2])
            
        self.likelihood = likelihood
        
    def sample(self, x, di):
        if self.nglobal_parameters > 0:
            return self.sampling_function(x, self.theta_g[di], self.boundaries, self.theta[di][self.partition_map])
        else:
            return self.sampling_function(x, self.boundaries, self.theta[di][self.partition_map])
    
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
    
    def _update_partition(self, datasets, di, pi, xi, xj):
        
        # Compute OLS and determine model probabilities
        prob = utils.update_partition_single_model(datasets[di],
                                                   self.gf_values[di],
                                                   self.lf_all_domain[di],
                                                   xi,
                                                   xj+1,
                                                   self.nlocal_parameters,
                                                   self._mu,
                                                   self._sigma,
                                                   self._S,
                                                   self.X[di],
                                                   self._N,
                                                   self._Q,
                                                   self._L,
                                                   self._U,
                                                   self._P,
                                                   self._aux1,
                                                   self._aux2,
                                                   self._auto_z,
                                                   self.prior_product[di],
                                                   self.B[di][pi],
                                                   self.Bi[di][pi],
                                                   self._u,
                                                   self.mu[di][pi],
                                                   self.theta[di][pi])
        
        self.ppratio[di][pi] = prob
        
        return prob
    
    def _update_partition_fixed(self, datasets, di, pi, xi, xj, ki, use_theta):
        
        prob = utils.compute_fixed_init(datasets[di],
                                        self.gf_values[di],
                                        self.lf_all_domain[di],
                                        xi,
                                        xj,
                                        self.kmax,
                                        self.nlocal_parameters,
                                        self._mu,
                                        self._sigma,
                                        self._S,
                                        self._detS,
                                        self._epsilon,
                                        self.X[di],
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
                                        self.pk[di][pi],
                                        self.kcdf[di][pi],
                                        self.prior_product[di][pi],
                                        self.B[di][pi],
                                        self.Bi[di][pi],
                                        self._u,
                                        self.k[di],
                                        pi,
                                        self.mu[di][pi],
                                        self.theta[di][pi],
                                        ki,
                                        use_theta)
        self.ppratio[di][pi] = prob
        
        return prob
