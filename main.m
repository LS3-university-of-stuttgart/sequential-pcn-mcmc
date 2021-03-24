%% Copyright 2020 Sebastian Reuschen    
% This is an example file to show the usage of the sequential pCN-MCMC. See
% the sequential_pcn_MCMC file for more information.
clear, clc, format compact, close all  

%% set flags: ONLY SET FLAGS IF YOU KNOW WHAT YOU ARE DOING. SETTING FLAGS INCORRECTLY WILL LEAD TO WRONG RESULTS!
flags = struct;
% flags.beta = 1;                                       % Sets beta to this value and does not optimize it internally.
% flags.kappa = 1;                                      % Sets kappa to this value and does not optimize it internally.
flags.parameter_tuning_frequency=250;                 % Frquency of tuning of the adaptive sequential pcn-MCMC

%% define prior model
model.discretization = [50,50];                         % discretization of unknowns
model.size = [5000,5000];                               % size of domain in meter
model.name = 'exponential';                             % geostatistical model for unknowns
model.variance   = 1;                                   % geostatistical parameter for exponential model
model.lambda     = [1500  2000];                        % length parameter in [x,y] direction of geostatistical model
model.rotation   = 135;                                 % Rotation of the lambda values (clovkwise) in degree.
model.mu         = -2.5;                                % Mean value of constant mean function.

model.id         = 1;                                   % id of this skript. For parallel execution, give each run a different id. 
model.saving_distance = 10;                             % Save theta and log_likelihood every saving_distance. Values larger one save Memory.
model.information_distance = 100;                       % Print log_likelihood, kappa, beta and acceptance rate to console every information_distance. 
model.visualization_distance = 100;                     % Plot current phi and simulation result every visualization_distance.
% Set any of them to 0 to not plot/show/save at all.


%% define log-likelihood function
solver.kriging = 1;                                                      % Use a Kriging example? If you want an example with a groundwater simulator. Please contact the authors of this code. 
solver.log_likelihood_function_handle = @compute_log_likelihood_example; % This defines the function which computes the log_likelihood. This function needs to be changed to suit your needs.
solver.set_up_solver_function_handle  = @set_up_solver_example;          % This function gets called ones before the start of the algorithm to set up the grid, reference solution or similar.


%% run MCMC
[flags,model,solver] = sequential_pcn_MCMC(50000,model,solver,flags);


%% save results
save('a_nice_name_to_save_a_file.mat','flags','model','solver','-v7.3');
disp('Results saved');