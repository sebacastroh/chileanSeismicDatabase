# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 13:48:02 2021

@author: sebac
"""
import utils
import numpy as np

class model1D:
    
    def __init__(self, npoints: np.ndarray, ndatasets: int, nglobal_parameters: int, kmax: int,
                 nlocal_parameters: np.ndarray, min_partitions: int, max_partitions: int,
                 xmin: float, xmax: float, pd: float, computeDesignMatrices: bool,
                 nprs=None, rs=None, modelType='heaviside') -> None:
        """
        

        Parameters
        ----------
        npoints : np.ndarray
            DESCRIPTION.
        ndatasets : int
            DESCRIPTION.
        nglobal_parameters : int
            DESCRIPTION.
        kmax : int
            DESCRIPTION.
        nlocal_parameters : np.ndarray
            DESCRIPTION.
        min_partitions : int
            DESCRIPTION.
        max_partitions : int
            DESCRIPTION.
        xmin : float
            DESCRIPTION.
        xmax : float
            DESCRIPTION.
        pd : float
            DESCRIPTION.
        computeDesignMatrices : bool
            DESCRIPTION.

        Returns
        -------
        None
            DESCRIPTION.

        """
        # Datasets parameters
        self.npoints = npoints.copy()
        self.ndatasets = ndatasets
        self.xmin = xmin
        self.xmax = xmax
        
        # Models parameters
        self.nglobal_parameters = nglobal_parameters
        self.kmax = kmax
        self.nlocal_parameters = nlocal_parameters.copy()
        
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
        self.global_parameters = np.empty((ndatasets, nglobal_parameters))
        self.nlmax = np.max(nlocal_parameters)
        
        # Global function values
        self.gf_values = [np.zeros(npoints[di]) for di in range(ndatasets)]
        
        # Design matrices
        self.computeDesignMatrices = computeDesignMatrices
        self.X = [np.zeros((kmax, npoints[di], self.nlmax)) for di in range(ndatasets)]
        
        # Local functions values
        self.lf_all_domain = [np.zeros(npoints[di]) for di in range(ndatasets)]
        self.lf_values = [np.zeros((npoints[di], self.max_partitions)) for di in range(ndatasets)]
        
        # Current local models
        self.k = [-np.ones(self.max_partitions, dtype=int) for di in range(ndatasets)] # Selected sub-models
        self.theta = [np.zeros((self.max_partitions, self.nlmax)) for di in range(ndatasets)] # Sampled values of local parameters
        self.mu = [np.zeros((self.max_partitions, self.nlmax)) for di in range(ndatasets)] # Mean values of local parameters
        self.B = [np.zeros((self.max_partitions, self.nlmax, self.nlmax)) for di in range(ndatasets)] # Cholesky decomposition of covariance matrices
        self.Bi = [np.zeros((self.max_partitions, self.nlmax, self.nlmax)) for di in range(ndatasets)] # Inverses of the Cholesky decomposition of covariance matrices
        self.pk = [np.zeros((self.max_partitions, kmax)) for di in range(ndatasets)] # Prior probability of sub-models using OLS
        self.kcdf = [np.zeros((self.max_partitions, kmax)) for di in range(ndatasets)] # Cumulative prior probabilities
        self.prior_product = [np.zeros((self.max_partitions, kmax)) for di in range(ndatasets)]
        self.ppratio = [np.zeros(self.max_partitions) for di in range(ndatasets)]
        
        # RJMCMC
        self.likelihood = 0.
        
        # For internal computations
        self._u = np.zeros(self.nlmax) # Random vector
        self._mu = np.zeros((kmax, self.nlmax)) # Mean of local parameters
        self._sigma = np.ones((kmax, self.nlmax)) # Standard deviation of local parameters
        self._S = np.zeros((kmax, self.nlmax, self.nlmax)) # Covariance matrices of each local model
        self._detS = np.zeros(kmax) # Determinants of covariance matrices
        self._B = np.zeros((self.nlmax, self.nlmax)) # Cholesky decomposition of covariance matrix
        self._Bi = np.zeros((self.nlmax, self.nlmax)) # Inverse of the Cholesky decomposition of covariance matrix
        self._epsilon = np.zeros(kmax) # Residual of each local model
        self._N = np.zeros((self.nlmax, self.nlmax)) # Normal matrix
        self._Q = np.zeros((self.nlmax, self.nlmax)) # Cofactor matrix (N^T*N)^-1
        self._L = np.zeros((self.nlmax, self.nlmax)) # Lower triangular matrix
        self._U = np.zeros((self.nlmax, self.nlmax)) # Upper triangular matrix
        self._P = np.zeros((self.nlmax, self.nlmax)) # Permutation matrix
        self._aux1 = np.zeros((self.nlmax, self.nlmax)) # Auxiliar matrix
        self._aux2 = np.zeros((self.nlmax, self.nlmax)) # Auxiliar matrix
        self._autoprior = np.zeros(kmax)
        self._alpha = np.zeros((kmax, kmax)) # Acceptance probability matrix
        self._auto_z = 3.
        
        self.nprs = nprs
        self.rs = rs
        self.modelType = modelType
    
    def initialize(self, datasets, global_parameters, global_function,
                   design_matrices, local_function):
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
                
            #Evaluate global function
            for di, dataset in enumerate(datasets):
                self.gf_values[di] = global_function(dataset[:,0], self.global_parameters[di])
        
        # Compute design matrices if corresponds
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                for ki in range(self.kmax):
                    self.X[di][ki,:,:self.nlocal_parameters[ki]] = design_matrices(dataset[:,0],
                                                                                   self.global_parameters[di],
                                                                                   self.boundaries,
                                                                                   ki,
                                                                                   self.nlocal_parameters[ki])
        else:
            for di in range(self.ndatasets):
                self.X[di] = design_matrices[di].copy()
        
        # Compute model for each segment
        for i, pi in enumerate(self.partition_map):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            
            for di, dataset in enumerate(datasets):
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                self._u = self.nprs.normal(size=self.nlmax)
                
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                self._update_partition(datasets, di, pi, xi, xj)
                
                ki = self.k[di][pi]
                nl = self.nlocal_parameters[ki]
                self.lf_values[di][:,pi] = local_function(dataset[:,0], self.global_parameters[di], xl, xr, ki, self.theta[di][pi,:nl])
        
    def initialize_fixed(self, datasets, fixed_init, global_parameters,
                         global_function, design_matrices, local_function):
        
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
                
            for di, dataset in enumerate(datasets):
                self.gf_values[di] = global_function(dataset[:,0], self.global_parameters[di])
            
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                for ki in range(self.kmax):
                    self.X[di][ki,:,:self.nlocal_parameters[ki]] = design_matrices(dataset[:,0],
                                                                                   self.global_parameters[di],
                                                                                   self.boundaries,
                                                                                   ki,
                                                                                   self.nlocal_parameters[ki])
        else:
            for di in range(self.ndatasets):
                self.X[di] = design_matrices[di].copy()
                
        # Compute model for each segment
        for i, pi in enumerate(self.partition_map):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            
            for di, dataset in enumerate(datasets):
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                self._u = self.nprs.normal(size=self.nlmax)
                
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                if 'models' in fixed_init and 'boundaries' in fixed_init:
                    ki = fixed_init['models'][di][pi]
                    nl = self.nlocal_parameters[ki]
                    if 'local_parameters' in fixed_init:
                        self.theta[di][pi,:nl] = fixed_init['local_parameters'][di][pi]
                        use_theta = 1
                    else:
                        use_theta = 0
                    
                    prob = self._update_partition_fixed(datasets, di, pi, xi, xj, ki, use_theta)
                else:
                    prob = self._update_partition(datasets, di, pi, xi, xj)
                
                if prob <= 0.:
                    return 0
                
                ki = self.k[di][pi]
                nl = self.nlocal_parameters[ki]
                self.lf_values[di][:,pi] = local_function(dataset[:,0], self.global_parameters[di], xl, xr, ki, self.theta[di][pi,:nl])
            
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
        self.nlmax = other.nlmax
        
        # Global function values
        self.gf_values = other.gf_values.copy()
        
        # Design matrices
        self.computeDesignMatrices = other.computeDesignMatrices
        self.X = [other.X[di].copy() for di in range(other.ndatasets)]
        
        # Local functions values
        self.lf_all_domain = [other.lf_all_domain[di].copy() for di in range(other.ndatasets)]
        self.lf_values = [other.lf_values[di].copy() for di in range(other.ndatasets)]
        
        # Current local models
        self.k = [other.k[di].copy() for di in range(other.ndatasets)]
        self.theta = [other.theta[di].copy() for di in range(other.ndatasets)]
        self.mu = [other.mu[di].copy() for di in range(other.ndatasets)]
        self.B = [other.B[di].copy() for di in range(other.ndatasets)]
        self.Bi = [other.Bi[di].copy() for di in range(other.ndatasets)]
        self.pk = [other.pk[di].copy() for di in range(other.ndatasets)]
        self.kcdf = [other.kcdf[di].copy() for di in range(other.ndatasets)]
        self.prior_product = [other.prior_product[di].copy() for di in range(other.ndatasets)]
        self.ppratio = [other.ppratio[di].copy() for di in range(other.ndatasets)]
        
        # RJMCMC
        self.likelihood = other.likelihood
        
    
    def perturb(self, process, datasets, design_matrices, local_function, global_function, global_parameters):
        
        if process == 0:
            status, process_prob = self.propose_birth(datasets, design_matrices, local_function)
        elif process == 1:
            status, process_prob = self.propose_death(datasets, design_matrices, local_function)
        elif process == 2:
            status, process_prob = self.propose_move(datasets, design_matrices, local_function)
        elif process == 3:
            status, process_prob = self.propose_merge(datasets, design_matrices, local_function)
        elif process == 4:
            status, process_prob = self.propose_local_parameters(datasets, local_function)
        else:
            status, process_prob = self.propose_global_parameters(datasets, global_parameters, global_function, design_matrices, local_function)
        
        return status, process_prob
    
    def propose_birth(self, datasets, design_matrices, local_function):
        
        if self.npartitions == self.max_partitions:
            return 0, 0.
        
        new_x = self.rs.random()*(self.xmax - self.xmin) + self.xmin
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
        
        if self.modelType == 'heaviside':
            for di in range(self.ndatasets):
                prob /= self.ppratio[di][self.partition_map[prev_i:]].prod()
        else:
            for di in range(self.ndatasets):
                prob /= self.ppratio[di][self.partition_map[prev_i]].prod()
        
        self.boundaries.insert(i, new_x)
        self.partition_map.insert(i, new_iy)
        self.npartitions += 1
        self.nboundaries += 1
        
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                for ki in range(self.kmax):
                    self.X[di][ki,:,:self.nlocal_parameters[ki]] = design_matrices(dataset[:,0],
                                                                                   self.global_parameters[di],
                                                                                   self.boundaries,
                                                                                   ki,
                                                                                   self.nlocal_parameters[ki])
        
        count = 0
        for i, pi in enumerate(self.partition_map[prev_i:], prev_i):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            for di, dataset in enumerate(datasets):
                
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                
                if pi == new_iy:
                    ppi = self.partition_map[prev_i]
                    ki = self.k[di][ppi]
                    nl = self.nlocal_parameters[ki]
                    self._u[:nl] = np.dot(self.Bi[di][ppi,:nl,:nl], self.theta[di][ppi,:nl] - self.mu[di][ppi,:nl])
                else:
                    ki = self.k[di][pi]
                    nl = self.nlocal_parameters[ki]
                    self._u[:nl] = np.dot(self.Bi[di][pi,:nl,:nl], self.theta[di][pi,:nl] - self.mu[di][pi,:nl])
                    
                if nl < self.nlmax:
                    self._u[nl:] = self.nprs.normal(size=self.nlmax - nl)
                    
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                curve_prob = self._update_partition(datasets, di, pi, xi, xj)
                if curve_prob <= 0:
                    return 0, 0.
                
                prob *= curve_prob
                
                ki = self.k[di][pi]
                nl = self.nlocal_parameters[ki]
                self.lf_values[di][:,pi] = local_function(dataset[:,0], self.global_parameters[di], xl, xr, ki, self.theta[di][pi,:nl])
            
            count += 1
            if self.modelType == 'pulse' and count == 2:
                break
        
        return 1, prob
                
    
    def propose_death(self, datasets, design_matrices, local_function):
        
        if self.npartitions <= self.min_partitions:
            return 0, 0.
        
        del_iy = self.nprs.choice(self.npartitions-1) + 1
        
        prev_i = del_iy-1
        
        prob = 1.
        
        if self.modelType == 'heaviside':
            for di in range(self.ndatasets):
                prob /= self.ppratio[di][self.partition_map[prev_i:]].prod()
        else:
            for di in range(self.ndatasets):
                prob /= self.ppratio[di][self.partition_map[prev_i:prev_i+2]].prod()
                
        self.partition_map.pop(del_iy)
        self.boundaries.pop(del_iy)
        self.npartitions -= 1
        self.nboundaries -= 1
        
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                for ki in range(self.kmax):
                    self.X[di][ki,:,:self.nlocal_parameters[ki]] = design_matrices(dataset[:,0],
                                                                                   self.global_parameters[di],
                                                                                   self.boundaries,
                                                                                   ki,
                                                                                   self.nlocal_parameters[ki])
        
        for i, pi in enumerate(self.partition_map[prev_i:], prev_i):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            for di, dataset in enumerate(datasets):
                
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                
                ki = self.k[di][pi]
                nl = self.nlocal_parameters[ki]
                self._u[:nl] = np.dot(self.Bi[di][pi,:nl,:nl], self.theta[di][pi,:nl] - self.mu[di][pi,:nl])
                    
                if nl < self.nlmax:
                    self._u[nl:] = self.nprs.normal(size=self.nlmax - nl)
                    
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                curve_prob = self._update_partition(datasets, di, pi, xi, xj)
                if curve_prob <= 0:
                    return 0, 0.
                
                prob *= curve_prob
                
                ki = self.k[di][pi]
                nl = self.nlocal_parameters[ki]
                self.lf_values[di][:,pi] = local_function(dataset[:,0], self.global_parameters[di], xl, xr, ki, self.theta[di][pi,:nl])
        
            if self.modelType == 'pulse':
                break
            
        return 1, prob
    
    def propose_move(self, datasets, design_matrices, local_function):
        
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
            if n <= self.nlmax:
                return 0, 0.
            n = len(np.where((dataset[:,0] >= new_x) & (dataset[:,0] <= xf))[0])
            if n <= self.nlmax:
                return 0, 0.
        
        prev_i = iy - 1
        self.boundaries[iy] = new_x
        
        prob = 1.
        
        if self.modelType == 'heaviside':
            for di in range(self.ndatasets):
                prob /= self.ppratio[di][self.partition_map[prev_i:]].prod()
        else:
            for di in range(self.ndatasets):
                prob /= self.ppratio[di][self.partition_map[prev_i:prev_i+2]].prod()
            
        if self.computeDesignMatrices:
            for di, dataset in enumerate(datasets):
                for ki in range(self.kmax):
                    self.X[di][ki,:,:self.nlocal_parameters[ki]] = design_matrices(dataset[:,0],
                                                                                   self.global_parameters[di],
                                                                                   self.boundaries,
                                                                                   ki,
                                                                                   self.nlocal_parameters[ki])
                    
        count = 0
        for i, pi in enumerate(self.partition_map[prev_i:], prev_i):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            for di, dataset in enumerate(datasets):
                
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                
                ki = self.k[di][pi]
                nl = self.nlocal_parameters[ki]
                self._u[:nl] = np.dot(self.Bi[di][pi,:nl,:nl], self.theta[di][pi,:nl] - self.mu[di][pi,:nl])
                    
                if nl < self.nlmax:
                    self._u[nl:] = self.nprs.normal(size=self.nlmax - nl)
                    
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                curve_prob = self._update_partition(datasets, di, pi, xi, xj)
                if curve_prob <= 0:
                    return 0, 0.
                
                prob *= curve_prob
                
                ki = self.k[di][pi]
                nl = self.nlocal_parameters[ki]
                self.lf_values[di][:,pi] = local_function(dataset[:,0], self.global_parameters[di], xl, xr, ki, self.theta[di][pi,:nl])
            
            count += 1
            if self.modelType == 'pulse' and count == 2:
                break
            
        return 1, prob
    
    def propose_merge(self, datasets, design_matrices, local_function):
        
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
        
        prob = 1.
        if self.modelType == 'heaviside':
            for di in range(self.ndatasets):
                prob /= self.ppratio[di][self.partition_map[prev_i:]].prod()
        else:
            for di in range(self.ndatasets):
                prob /= self.ppratio[di][self.partition_map[prev_i:prev_i+3]].prod()
        
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
                for ki in range(self.kmax):
                    self.X[di][ki,:,:self.nlocal_parameters[ki]] = design_matrices(dataset[:,0],
                                                                                   self.global_parameters[di],
                                                                                   self.boundaries,
                                                                                   ki,
                                                                                   self.nlocal_parameters[ki])
        
        count = 0
        for i, pi in enumerate(self.partition_map[prev_i:], prev_i):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            for di, dataset in enumerate(datasets):
                dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
                xi, xj = dataset_indices[[0, -1]]
                
                ki = self.k[di][pi]
                nl = self.nlocal_parameters[ki]
                self._u[:nl] = np.dot(self.Bi[di][pi,:nl,:nl], self.theta[di][pi,:nl] - self.mu[di][pi,:nl])
                    
                if nl < self.nlmax:
                    self._u[nl:] = self.nprs.normal(size=self.nlmax - nl)
                    
                self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
                
                curve_prob = self._update_partition(datasets, di, pi, xi, xj)
                if curve_prob <= 0:
                    return 0, 0.
                
                prob *= curve_prob
                
                ki = self.k[di][pi]
                nl = self.nlocal_parameters[ki]
                self.lf_values[di][:,pi] = local_function(dataset[:,0], self.global_parameters[di], xl, xr, ki, self.theta[di][pi,:nl])
        
            count += 1
            if self.modelType == 'pulse' and count == 3:
                break
        
        return 1, prob
    
    def propose_local_parameters(self, datasets, local_function):
        
        if self.ndatasets == 1:
            di = 0
        else:
            di = self.nprs.choice(self.ndatasets)
            
        iy = self.nprs.choice(self.npartitions)
        
        if self.modelType == 'heaviside':
            prob = 1./self.ppratio[di][self.partition_map[iy:]].prod()
        else:
            prob = 1./self.ppratio[di][self.partition_map[iy]]
        
        for i, pi in enumerate(self.partition_map[iy:], iy):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            dataset = datasets[di]
            dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
            xi, xj = dataset_indices[[0, -1]]
            
            ki = self.k[di][pi]
            nl = self.nlocal_parameters[ki]
            self._u[:nl] = np.dot(self.Bi[di][pi,:nl,:nl], self.theta[di][pi,:nl] - self.mu[di][pi,:nl])
                
            if nl < self.nlmax:
                self._u[nl:] = self.nprs.normal(size=self.nlmax - nl)
                
            self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
            
            curve_prob = self._update_partition(datasets, di, pi, xi, xj)
            if curve_prob <= 0:
                return 0, 0.
            
            prob *= curve_prob
            
            ki = self.k[di][pi]
            nl = self.nlocal_parameters[ki]
            self.lf_values[di][:,pi] = local_function(dataset[:,0], self.global_parameters[di], xl, xr, ki, self.theta[di][pi,:nl])
            
            if self.modelType == 'pulse':
                break
        
        return 1, prob
    
    def propose_global_parameters(self, datasets, global_parameters, global_function, design_matrices, local_function):
        
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
        
        self.gf_values[di] = global_function(datasets[di][:,0], self.global_parameters[di])
        
        if self.computeDesignMatrices:
            for ki in range(self.kmax):
                self.X[di][ki,:,:self.nlocal_parameters[ki]] = design_matrices(datasets[di][:,0],
                                                                               self.global_parameters[di],
                                                                               self.boundaries,
                                                                               ki,
                                                                               self.nlocal_parameters[ki])
        
        prob = 1./self.ppratio[di][self.partition_map].prod()
        
        for i, pi in enumerate(self.partition_map):
            xl = self.boundaries[i]
            xr = self.boundaries[i+1]
            dataset = datasets[di]
            dataset_indices = np.where((dataset[:,0] >= xl) & (dataset[:,0] <= xr))[0]
            xi, xj = dataset_indices[[0, -1]]
            
            ki = self.k[di][pi]
            nl = self.nlocal_parameters[ki]
            self._u[:nl] = np.dot(self.Bi[di][pi,:nl,:nl], self.theta[di][pi,:nl] - self.mu[di][pi,:nl])
                
            if nl < self.nlmax:
                self._u[nl:] = self.nprs.normal(size=self.nlmax - nl)
                
            self.lf_all_domain[di] = self.lf_values[di][:,self.partition_map[:i]].sum(1)
            
            curve_prob = self._update_partition(datasets, di, pi, xi, xj)
            if curve_prob <= 0:
                return 0, 0.
            
            prob *= curve_prob
            
            ki = self.k[di][pi]
            nl = self.nlocal_parameters[ki]
            self.lf_values[di][:,pi] = local_function(dataset[:,0], self.global_parameters[di], xl, xr, ki, self.theta[di][pi,:nl])
        
        return 1, prob
        
        
    
    def misfit(self, datasets: list, likelihood_function):
        
        likelihood = 0.
        for di, dataset in enumerate(datasets):
            likelihood += likelihood_function(self.global_parameters[di],
                                              self.boundaries,
                                              self.k[di][self.partition_map],
                                              self.theta[di][self.partition_map],
                                              self.nlocal_parameters,
                                              dataset[:,0],
                                              dataset[:,1],
                                              dataset[:,2])
            
        self.likelihood = likelihood
    
    def _check_data(self, i, x, datasets):
        xi = self.boundaries[i-1]
        xf = self.boundaries[i]
        
        for dataset in datasets:
            n = len(np.where((dataset[:,0] >= xi) & (dataset[:,0] < x))[0])
            if n <= self.nlmax:
                return False
            n = len(np.where((dataset[:,0] >= x) & (dataset[:,0] <= xf))[0])
            if n <= self.nlmax:
                return False
        
        return True
    
    def _update_partition(self, datasets, di, pi, xi, xj):
        
        # Compute OLS and determine model probabilities
        prob = utils.update_partition(datasets[di],
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
    