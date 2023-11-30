# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 10:43:45 2022

@author: sebac
"""

import numpy as np

class mp_so_connected_lines:

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
        
        self.ymin = np.empty(ndatasets).tolist()
        self.ymax = np.empty(ndatasets).tolist()
        self.dy = np.empty(ndatasets).tolist()
        
        self.x = [np.empty(npo) for npo in self.npoints]
        self.y = [np.empty(npo) for npo in self.npoints]
        self.s = [np.empty(npo) for npo in self.npoints]
        
        self.initiated = False

        # Models parameters
        self.nglobal_parameters = 0 # nglobal_parameters
        self.nlocal_parameters = 2 # nlocal_parameters

        # Partition parameters
        self.npartitions = 0
        self.min_partitions = min_partitions
        self.max_partitions = max_partitions
        self.partition_map = [0]
        
        # For the sake of simplicity, we define boundaries
        self.nboundaries = 0
        self.boundaries = []
        self.yboundaries = [[] for di in range(ndatasets)]
        self.locations = [[] for di in range(ndatasets)]
        self.min_boundaries = min_partitions + 1
        self.max_boundaries = max_partitions + 1
        
        # Boundary perturbation
        self.pd = pd

        # Current local models
        self.theta = [np.empty((self.max_partitions, self.nlocal_parameters)) for di in range(ndatasets)] # Sampled values of local parameters
        
        self.nprs = nprs

        # RJMCMC
        self.likelihood = 0.

        # Functions
        self.sampling_function = sampling_function
        self.likelihood_function = loglikelihood_function
        self.global_function = global_function
        self.local_function = local_function
        
        self.samples = samples

    def initialize(self, datasets):
        # By default the first model uses the minimum partitions allowed
        self.initiated = True
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

        # Compute model for each segment
        boundaries = np.array(self.boundaries)
        difference = (boundaries[1:] - boundaries[:-1])
        for di, dataset in enumerate(datasets):
            self.ymin[di] = dataset[:,1].min()
            self.ymax[di] = dataset[:,1].max()
            self.dy[di] = self.ymax[di] - self.ymin[di]
            
            self.x[di] = np.ascontiguousarray(dataset[:,0])
            self.y[di] = np.ascontiguousarray(dataset[:,1])
            self.s[di] = np.ascontiguousarray(dataset[:,2])
            
            sampledYBoundaries = self.nprs.uniform(self.ymin[di], self.ymax[di], self.nboundaries)
            self.yboundaries[di] = sampledYBoundaries.tolist()
            
            self.locations[di] = []
            for x in self.boundaries[:-1]:
                pos = np.argmax(self.x[di] >= x)
                self.locations[di].append(pos)
            self.locations[di].append(self.npoints[di])
            
            a = (sampledYBoundaries[1:] - sampledYBoundaries[:-1])/difference
            self.theta[di][self.partition_map,0] = a
            self.theta[di][self.partition_map,1] = sampledYBoundaries[1:] - a*boundaries[1:]
            

    def initialize_fixed(self, datasets, fixed_init):
        pass

    def clone(self, other):
        
        if not self.initiated:
            self.ymin = other.ymin.copy()
            self.ymax = other.ymax.copy()
            self.dy = other.dy.copy()
            self.x = [other.x[di].copy() for di in range(self.ndatasets)]
            self.y = [other.y[di].copy() for di in range(self.ndatasets)]
            self.s = [other.s[di].copy() for di in range(self.ndatasets)]
            
            self.initiated = True
        
        # Partition parameters
        self.npartitions = other.npartitions
        self.min_partitions = other.min_partitions
        self.max_partitions = other.max_partitions
        self.partition_map = other.partition_map.copy()
        
        # For the sake of simplicity, we define boundaries
        self.nboundaries = other.nboundaries
        self.boundaries = other.boundaries.copy()
        self.yboundaries = [ybounds.copy() for ybounds in other.yboundaries]
        self.locations = [locations.copy() for locations in other.locations]
        self.min_boundaries = other.min_boundaries
        self.max_boundaries = other.max_boundaries
        
        # Boundary perturbation
        self.pd = other.pd
        
        self.nglobal_parameters = 0        
        self.nlocal_parameters = 2
        
        # Current local models
        self.theta = [theta.copy() for theta in other.theta]
        
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
        else:
            status, process_prob = self._propose_local_parameters(datasets)
            
        boundaries = np.array(self.boundaries)
        difference = (boundaries[1:] - boundaries[:-1])
        
        for di, dataset in enumerate(datasets):
            sampledYBoundaries = np.array(self.yboundaries[di])
            a = (sampledYBoundaries[1:] - sampledYBoundaries[:-1])/difference
            b = sampledYBoundaries[1:] - a*boundaries[1:]
            self.theta[di][self.partition_map,0] = a
            self.theta[di][self.partition_map,1] = b
        
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
        
        pos = self._check_data(i, new_x, datasets)
        if not pos:
            return 0, 0.
        
        self.boundaries.insert(i, new_x)
        self.partition_map.insert(i, new_iy)
        self.npartitions += 1
        self.nboundaries += 1
        
        prob = 1.
        
        for di, dataset in enumerate(datasets):
            sampledYBoundary = self.nprs.uniform(self.ymin[di], self.ymax[di])
            self.yboundaries[di].insert(i, sampledYBoundary)
            self.locations[di].insert(i, pos[di])
            prob /= self.dy[di]
        
        return 1, prob
    
    def _propose_death(self, datasets):
        
        if self.npartitions <= self.min_partitions:
            return 0, 0.
        
        del_iy = self.nprs.integers(0, self.npartitions-1) + 1
                
        self.partition_map.pop(del_iy)
        self.boundaries.pop(del_iy)
        
        prob = 1.
        
        for di in range(self.ndatasets):
            self.yboundaries[di].pop(del_iy)
            self.locations[di].pop(del_iy)
            prob *= self.dy[di]
            
        self.npartitions -= 1
        self.nboundaries -= 1
            
        return 1, prob
    
    def _propose_move(self, datasets):
        
        if self.npartitions <= self.min_partitions:
            return 0, 0.
        
        iy = self.nprs.integers(0, self.npartitions-1) + 1
        
        old_x = self.boundaries[iy]
        new_x = old_x + self.nprs.normal(0., self.pd)
        
        xi = self.boundaries[iy-1]
        xf = self.boundaries[iy+1]
        
        if new_x <= xi or new_x >= xf:
            return 0, 0.
        
        pos = self._check_data(iy, new_x, datasets)
        if not pos:
            return 0, 0.
        
        self.boundaries[iy] = new_x
        for di, dataset in enumerate(datasets):
            self.locations[di][iy] = pos[di]
        
        prob = 1.
            
        return 1, prob
    
    def _propose_merge(self, datasets):
        
        if self.npartitions <= self.min_partitions or self.nboundaries <= 3:
            return 0, 0.
        
        jy = self.nprs.integers(0, self.npartitions-1) + 1
        
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
        ya = [ybounds[jy] for ybounds in self.yboundaries]
        if case == 1:
            ik = 1
            xb = self.boundaries[ky]
            yb = [ybounds[ky] for ybounds in self.yboundaries]
        elif case == 2:
            ik = 0
            xb = self.boundaries[iy]
            yb = [ybounds[iy] for ybounds in self.yboundaries]
        else:
            ik = self.nprs.integers(0, 2)
            if ik == 0:
                xb = self.boundaries[iy]
                yb = [ybounds[iy] for ybounds in self.yboundaries]
            else:
                xb = self.boundaries[ky]
                yb = [ybounds[ky] for ybounds in self.yboundaries]
        
        alpha = self.nprs.random()
        new_x = xa*alpha + xb*(1.-alpha)
        
        if ik == 0:
            self.boundaries[iy] = new_x
        else:
            self.boundaries[ky] = new_x
            
        self.partition_map.pop(jy)
        self.boundaries.pop(jy)
        
        prob = 1.
            
        for di, dataset in enumerate(datasets):
            new_y = ya[di]*alpha + yb[di]*(1.-alpha)
            p1 = self.locations[di][jy-1]
            p2 = self.locations[di][jy+1]
            pos = np.searchsorted(self.x[di][p1:p2], new_x) + p1
            if ik == 0:
                self.yboundaries[di][iy] = new_y
                self.locations[di][iy] = pos
            else:
                self.yboundaries[di][ky] = new_y
                self.locations[di][ky] = pos
                
            self.yboundaries[di].pop(jy)
            self.locations[di].pop(jy)
            
            prob *= self.dy[di]
            
            
        self.npartitions -= 1
        self.nboundaries -= 1
            
        return 1, prob
    
    def _propose_local_parameters(self, datasets):
        
        if self.ndatasets == 1:
            di = 0
        else:
            di = self.nprs.integers(0, self.ndatasets)
        
        bi = self.nprs.integers(0, self.nboundaries)
        sampledYBoundary = self.nprs.uniform(self.ymin[di], self.ymax[di])
        self.yboundaries[di][bi] = sampledYBoundary
        
        prob = 1.
        
        return 1, prob
    
    def _check_data(self, i, x, datasets):
        pos = []
        
        for di, dataset in enumerate(datasets):
            p1 = self.locations[di][i-1]
            p2 = self.locations[di][i]
            if p1 == p2:
                return False
            
            p = np.searchsorted(self.x[di][p1:p2], x) + p1
            
            if (p - p1) <= 2:
                return False
            elif (p2 - p) <= 2:
                return False
            
            pos.append(p)
        
        return pos
    
    def _sampling_function(self, x, locations, local_values):
        y = np.empty_like(x)
        
        for i in range(len(locations)-1):
            p1 = locations[i]
            p2 = locations[i+1]
            y[p1:p2] = local_values[i,0]*x[p1:p2] + local_values[i,1]
        
        return y
    
    def _likelihood_function(self, locations, local_values, x, data, sigma):
        y = self._sampling_function(x, locations, local_values)
        phi = np.sum(((y - data)/(2.*sigma))**2)
        
        return phi
    
    def misfit(self, datasets):
        likelihood = 0.
        for di, dataset in enumerate(datasets):
            likelihood += self._likelihood_function(self.locations[di],
                                                    self.theta[di][self.partition_map],
                                                    self.x[di],
                                                    self.y[di],
                                                    self.s[di])
            
        self.likelihood = likelihood

    def sample(self, x, di):
        locations = []
        for xi in self.boundaries[:-1]:
            pos = np.argmax(x >= xi)
            locations.append(pos)
        locations.append(self.samples)
        
        return self._sampling_function(x, locations, self.theta[di][self.partition_map])
