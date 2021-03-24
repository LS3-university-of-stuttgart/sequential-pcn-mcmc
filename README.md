# README #

This repository implements the [sequential pCN-MCMC](https://bitbucket.org/Reuschen/sequential-pcn-mcmc) as presented in the reference below. 

### Usage ###

The sequential pCN-MCMC is implemented in the *sequenial_pcn_MCMC.m* file. This package comes with a Kriging example of the sequential pCN-MCMC. By Kriging example, we understand a case where we infer a spatially distributed variable by measuring the variable directly. Here, we measure and infer the hydraulic conductivity.

More examples, as described in the publication are available by contacting the authors. In these examples, the hydraulic conductivity is inferred by measuring the hydraulic head. Thereby, we use a FEM solver to calculate the hydraulic heads based on the hydraulic conductivity. This FEM solver is not open source and we cannot provide it here.

To run the artificial Kriging example, run the *main.m* file. This code is implemented and tested in MATLAB2018a and MATLAB 2020a. It may not work properly in other version.

### Setting up your model ###
 
If you want to combine the sequential pCN-MCMC with your own model, you need to adapt the *compute_log_likelihood_example.m* and *set_up_solver_example.m* accordingly. The *set_up_solver_example.m* file is called once before the optimization. This file can be used to set up your simulation code (e.g. set up the Mesh, set boundary conditions, ...). The *compute_log_likelihood_example.m* is called every time the likelihood is evaluated. Make sure that this file is as efficient as possible (e.g. do not set up the mesh here). Of course, you can change the names of these two files as long as you change them in the *main.m* file as well.

You can define the variogram (covariance matrix) type and parameters of your prior model in the *main.m* file. Currently, exponential and Matern (with &nu; = 2.5) covariance’s are available. You can add more covariance models in the *sequential_pcn_MCMC.m* file in the *compute_covariance_matrix* function. 

### Reference ###
If you use this algorithm for your research, please reference the article  

> Sebastian Reuschen, Fabian Jobst and Wolfgang Nowak  
> **Sequential pCN-MCMC, an efficient MCMC method for Bayesian inversion of high-dimensional multi-Gaussian priors**  
> *Water Resources Research* 

### Copyright ###

Copyright 2020 Sebastian Reuschen                                         

See the file COPYING and COPYING.LESSER for full copying permissions.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.                                    
                                                                      
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.                           
                                                                       
You should have received a copy of the GNU Lesser General Public License along with this program. If not, see <http://www.gnu.org/licenses/>. 

### Contact ###
For more information contact sebastian.reuschen@iws.uni-stuttgart.de or wolfgang.nowak@iws.uni-stuttgart.de
