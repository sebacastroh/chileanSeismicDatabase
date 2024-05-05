# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 11:38:10 2021

@author: sebac
"""
import numpy as np

class result1D:
    """
    """
    def __init__(self, ndatasets, kmax, max_partitions,
                 burnin, total, nglobal_parameters, nlocal_parameters,
                 samples, xmin, xmax):
        """
        

        Parameters
        ----------
        ndatasets : int
            DESCRIPTION.
        burnin : int
            DESCRIPTION.
        total : int
            DESCRIPTION.
        nglobal_parameters : int
            DESCRIPTION.
        nlocal_parameters : np.ndarray
            DESCRIPTION.
        samples : int
            DESCRIPTION.
        xmin : float
            DESCRIPTION.
        xmax : float
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # Save inputs
        self.ndatasets = ndatasets
        self.kmax = kmax
        self.max_partitions = max_partitions
        
        self.burnin = burnin
        self.total = total
        
        self.nglobal_parameters = nglobal_parameters
        
        if isinstance(nlocal_parameters, int):
            self.nlocal_parameters = nlocal_parameters
        else:
            self.nlocal_parameters = nlocal_parameters.copy()
        
        self.samples = samples
        
        # Create arrays
        self.global_parameters = []
        self.models = []
        self.local_parameters = []
        
        self.misfit = []
        
        if self.max_partitions > 1:
            self.npartitions = []
            self.boundaries = []
            self.boundaries_histogram = None
            self.partitions_histogram = None
        
        self.acceptance_history = []
        self.processes_history = []
        
        self.y = []
        self.x = np.linspace(xmin, xmax, samples)
        
    def _sample(self, model, y, process, acceptance):
        
        self.misfit.append(model.likelihood)
        
        if self.max_partitions > 1:
            self.npartitions.append(model.npartitions)
            self.boundaries.append(model.boundaries)
        
        # Save global parameters
        if self.nglobal_parameters > 0:
            self.global_parameters.append(model.theta_g)
        
        # Save k
        if self.kmax > 1:
            if self.max_partitions == 1:
                current_k = model.k[0]
            else:
                if self.ndatasets == 1:
                    current_k = model.k.copy()
                else:
                    current_k = []
                    for di in range(self.ndatasets):
                        current_k.append(model.k[di][model.partition_map].tolist())
        
            self.models.append(current_k)
        
        # Save local parameters
        if self.max_partitions == 1:
            self.local_parameters.append(model.theta)
        else:
            if self.ndatasets == 1:
                self.local_parameters.append(model.theta[0][model.partition_map])
            else:    
                current_theta = []
                for di in range(self.ndatasets):
                    
                    if self.kmax == 1:
                        this_theta = model.theta[di][model.partition_map]
                    else:
                        this_theta = model.theta[di][model.partition_map].tolist()
                        nls = model.nlocal_parameters[model.k[di][model.partition_map]]
                        for i in range(len(this_theta)):
                            this_theta[i] = this_theta[i][:nls[i]]
                            
                    current_theta.append(this_theta)
                
                if self.kmax == 1:
                    current_theta = np.array(current_theta)
                    
                self.local_parameters.append(current_theta)
        
        if y is not None:
            if self.ndatasets == 1:
                self.y.append(y[0])
            else:
                self.y.append([this_y.copy() for this_y in y])
            
        self.processes_history.append(process)
        self.acceptance_history.append(acceptance)
            
    def _assemble(self, nprocesses):
        
        self.misfit = np.array(self.misfit)
        
        if self.max_partitions == 1:
            self.models = np.array(self.models)
        
        if self.max_partitions > 1:
            self.npartitions = np.array(self.npartitions)
        
        if self.nglobal_parameters > 0:
            self.global_parameters = np.array(self.global_parameters)
        
        self.y = np.array(self.y)
        
        if self.max_partitions == 1 and self.kmax == 1:
            self.local_parameters = np.array(self.local_parameters)
        
        if self.max_partitions > 1:
            bound_hist = []
            for boundaries in self.boundaries[self.burnin:]:
                bound_hist.extend(boundaries[1:-1])
                
            bound_hist, bound_edges = np.histogram(bound_hist,
                                                   bins=self.samples,
                                                   range=(self.x[0], self.x[-1]))
            
            self.boundaries_histogram = ((bound_edges[1:] + bound_edges[:-1])/2., bound_hist)
            
            partitions, count = np.unique(self.npartitions[self.burnin:], return_counts=True)
            self.partitions_histogram = (partitions, count)
        
        self.processes_history = np.array(self.processes_history)
        self.acceptance_history = np.array(self.acceptance_history)
        
        self.total_accepted = self.acceptance_history.sum()
        self.total_rejected = self.total - self.total_accepted
        
        if nprocesses > 1:
            self.processes = []
            self.acceptance_by_process = []
            self.rejection_by_process = []
            
            for i in range(nprocesses):
                pos = np.where(self.processes_history == i)[0]
                n = len(pos)
                self.processes.append(n)
                
                m = len(np.where(self.acceptance_history[pos] == 1)[0])
                self.acceptance_by_process.append(m)
                self.rejection_by_process.append(n-m)
                
            self.processes = np.array(self.processes)
            self.acceptance_by_process = np.array(self.acceptance_by_process)
            self.rejection_by_process = np.array(self.rejection_by_process)
