function [flags,model,solver] = sequential_pcn_MCMC(number_of_samples,model,solver,flags)  
%    The sequential_pcn_MCMC samples the posterior distribution defined by
%    the prior distribution in "model" and the the lieklihood distribution
%    implied by and evaluated by the "solver".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright 2020 Sebastian Reuschen                                         
% 
%    See the file COPYING and COPYING.LESSER for full copying permissions.
%
%    This program is free software: you can redistribute it and/or modify   
%    it under the terms of the GNU Lesser General Public License as published by   
%    the Free Software Foundation, either version 3 of the License, or      
%    (at your option) any later version.                                    
%                                                                          
%    This program is distributed in the hope that it will be useful,        
%    but WITHOUT ANY WARRANTY; without even the implied warranty of         
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           
%    GNU Lesser General Public License for more details.                           
%                                                                           
%    You should have received a copy of the GNU Lesser General Public License      
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This code was published with the paper "Reuschen et al. (2020). 
%   Sequential-pCN-MCMC, an efficient MCMC for Bayesian inversion of 
%   high-dimensional multi-Gaussian priors. Water Resources Research TODO
%   add number"
%   Please always cite this paper together with this code. 


    if(~exist('flags','var'))
        flags = struct;
    end
    randomize(model)
    [model,solver,flags] = check_inputs(number_of_samples,model,solver,flags);

    for i=1:number_of_samples
        %% Core algorithm
        if( flags.number_of_samples_max<flags.current_sample_number)
            disp('Maximum number of samples reached.');
            flags.finished = true;
            break
        end

        % Proposal
        flags.current_sample_number = flags.current_sample_number+1;
        [idx_1,idx_2]  = get_box_algorithm(flags,model);
        [mu_tilde,covariance_tilde] = compute_mean_and_covariance(idx_1,idx_2,flags);

        theta_tilde = flags.theta;
        if(numel(mu_tilde>0)) %It can happen, that the box is so small that nothing changes. This if clause catches this case.
            zeta = mvnrnd(zeros(numel(idx_1),1),covariance_tilde)';
            theta_tilde(idx_1) = sqrt(1-flags.beta.^2)* (flags.theta(idx_1) - mu_tilde  ) + flags.beta*zeta + mu_tilde ;
        end

        % Acceptance or Rejection step
        log_likelihood_tilde = solver.log_likelihood_function_handle(theta_tilde,model,solver,flags);
        alpha = min(exp(log_likelihood_tilde-flags.log_likelihood),1);
        r = rand(1);
        if r< alpha  %Accept
            flags.log_likelihood = log_likelihood_tilde;
            flags.theta = theta_tilde;
            flags.number_accepted = flags.number_accepted +1;
        end

        %% Saving and plotting
        if(mod(flags.current_sample_number,model.saving_distance)==0)
            flags.theta_all(:,flags.current_sample_number/model.saving_distance) = flags.theta;
            flags.log_likelihood_all(:,flags.current_sample_number/model.saving_distance) = flags.log_likelihood;
        end
        if(mod(flags.current_sample_number,model.information_distance)==0)
            disp([datestr(now,'YYYY_mm_DD_HH:MM:SS') ' id' num2str(model.id) ' Sample numer ' num2str(flags.current_sample_number)  ': log_likelihood = ' num2str(round(flags.log_likelihood,4,'significant')) ', acceptance_rate = '  num2str(round(flags.number_accepted/flags.current_sample_number,4,'significant')*100) '%, beta = ' num2str(flags.beta) ', kappa = ' num2str(flags.kappa)]);
        end
        if(mod(flags.current_sample_number,model.visualization_distance)==0)
            figure(2)
            subplot(1,2,1)        
            imagesc(model.x_pos,model.y_pos,reshape(flags.theta,[model.discretization(2),model.discretization(1)]));caxis([-5,0]);set(gca,'YDir','normal');daspect([1 1 1]);colorbar;
            xlabel('x [m]'); ylabel('y [m]'); title('theta');
            subplot(1,2,2); solver.log_likelihood_function_handle(flags.theta,model,solver,flags,1); %Inefficient because the function gets called again here.
            try
                suptitle(['Sample ' num2str(flags.current_sample_number)]);
            catch
                try
                    sgtitle(['Sample ' num2str(flags.current_sample_number)]);
                catch
                end
            end
            drawnow;
        end

       %% Apaptive parameter tuning
        if(flags.enable_parameter_tuning)
            flags.theta_recently(:,mod(flags.current_sample_number-1,flags.parameter_tuning_frequency)+1) = flags.theta;
            if(mod(flags.current_sample_number,flags.parameter_tuning_frequency)==0)
                flags = parameter_tuning(flags,model);
            end
        end

    end

end

%% Additional functions for the sequential pCN-MCMC.
function flags = parameter_tuning(flags,model)

    % Evaluate 
    index = flags.current_sample_number/flags.parameter_tuning_frequency; 
    flags.objective_function(ceil(index/4),mod(index-1,4)+1) = evaluate_hyper_parameters(flags,model);

    
    if(mod(index-1,4)+1==4) %after four different hyper-parameters were tested, choose new hyper-parameters.
        
        % Evaluate derivatives of beta and kappa
        diff_beta = flags.objective_function(ceil(index/4),1)- flags.objective_function(ceil(index/4),2);
        diff_kappa = flags.objective_function(ceil(index/4),3)- flags.objective_function(ceil(index/4),4);
        
        if(diff_beta == 0 &&diff_kappa==0) %This happens if nothing gets accepted or (very unlikely) by chance.
            diff_beta = -1;diff_kappa=-1;  %Walk to the bootom left. Otherwise the numerics will not work.
        end
        
        %If beta or kappa is user defined then set derivative to zero.
        if(flags.beta_is_constant)
            diff_beta = 0;
        end
        if(flags.kappa_is_constant)
            diff_kappa = 0;
        end

        %Compute how far to jump
        opt_count = ceil(index/4);
        if(ceil(index/4)>5) %10 %TODO
            beta_derivation = flags.all_tested_beta(1:opt_count-2) ./ flags.all_tested_beta(opt_count-1);
            kappa_derivation = flags.all_tested_kappa(1:opt_count-2) ./ flags.all_tested_kappa(opt_count-1);
         
            kappa_smaller = flags.delta_kappa_proportion^1.25>max(kappa_derivation,kappa_derivation.^-1);
            beta_smaller = flags.delta_beta_proportion^1.25>max(beta_derivation,beta_derivation.^-1);
            if( sum(kappa_smaller .* beta_smaller)>4)
                flags.jump_width_factor = 0.5;
            else
                flags.jump_width_factor = 1;
            end            
        else 
            flags.jump_width_factor = 2;
        end
        
        %Compute new beta and kappa positions
        flags.tested_beta  = min(max(flags.tested_beta   * flags.delta_beta_proportion  .^  (flags.jump_width_factor*diff_beta /sqrt(diff_beta^2+diff_kappa^2))  ,0),1);
        flags.tested_kappa  = min(max(flags.tested_kappa   * flags.delta_kappa_proportion  .^  (flags.jump_width_factor*diff_kappa /sqrt(diff_beta^2+diff_kappa^2))  ,1/(2*max(model.discretization))),1);
      
        flags.all_tested_beta(ceil(index/4))  = flags.tested_beta;
        flags.all_tested_kappa(ceil(index/4)) = flags.tested_kappa;

        if(flags.debug)
            objective_function_plot = flags.objective_function(ceil(index/4),:)/max(flags.objective_function(ceil(index/4),:))*100;
            disp(['^ kappa']);
            disp(['|          ' num2str(round(objective_function_plot(3),2))]);
            disp(['|' num2str(round(objective_function_plot(2),2)) '                 ' num2str(round(objective_function_plot(1),2)) ]);
            disp(['|          ' num2str(round(objective_function_plot(4),2))]);
            disp(['----------> beta']);
        end
        disp([ 'New values: id' num2str(model.id) ' tested beta = ' num2str(flags.tested_beta) ' , tested kappa = ' num2str(flags.tested_kappa) ,' , jump_width_factor = ', num2str(flags.jump_width_factor)]);
         
        % Test if stop tuning criterion is fulfilled
        opt_count = ceil(index/4);
        if(opt_count>1)
            beta_derivation = flags.all_tested_beta(1:opt_count-1) ./ flags.all_tested_beta(opt_count);
            kappa_derivation = flags.all_tested_kappa(1:opt_count-1) ./ flags.all_tested_kappa(opt_count);
            kappa_finish = flags.delta_kappa_proportion^0.75>max(kappa_derivation,kappa_derivation.^-1);
            beta_finish = flags.delta_beta_proportion^0.75>max(beta_derivation,beta_derivation.^-1);
            clear flags.theta_recently
            if(sum(kappa_finish .* beta_finish)>6|| (ceil(index/4)==flags.max_tuning_iterations)) 
                flags.kappa = flags.tested_kappa;
                flags.beta = flags.tested_beta;
                flags.enable_parameter_tuning = false;
                flags.end_of_parameter_tuning = flags.current_sample_number;
                disp('Finished parameter tuning.');
                if((mod(index-1,4)+1==4&&(ceil(index/4)==flags.max_tuning_iterations)))
                    disp('Care!: Tuning ended because maximum number of tuning iterations was reached! The optimal hyperparamter might not have been found.');
                end
                return
            end
        end
    end
    
    % Choose new beta and kappa
    if(mod(index-1,4)+1==4)
        flags.kappa = flags.tested_kappa;
        flags.beta  = min(max(flags.tested_beta  * flags.delta_beta_proportion ,0),1);  
    elseif(mod(index-1,4)+1==1)
        flags.kappa = flags.tested_kappa;
        flags.beta  = min(max(flags.tested_beta  * 1/flags.delta_beta_proportion ,0),1);    
    elseif(mod(index-1,4)+1==2)
        flags.beta = flags.tested_beta;
        flags.kappa  = min(max(flags.tested_kappa  * flags.delta_kappa_proportion ,1/(2*max(model.discretization))),1); 
    elseif(mod(index-1,4)+1==3)
        flags.beta = flags.tested_beta;
        flags.kappa  = min(max(flags.tested_kappa  * 1/flags.delta_kappa_proportion ,1/(2*max(model.discretization))),1);
    end
end

function    covariance = compute_covariance_matrix(pos_1,pos_2,model)
    if (strcmp(model.name,'exponential'))
        rotation_matrix = [cosd(model.rotation) -sind(model.rotation); sind(model.rotation) cosd(model.rotation)];
        lambda_squared = [(model.lambda(1).^2);(model.lambda(2).^2)];
        inside = sum((rotation_matrix*(pos_1'-pos_2')).^2./lambda_squared).^0.5;
        covariance = exp((-inside)).*model.variance; 
    elseif (strcmp(model.name,'matern2.5'))
        rotation_matrix = [cosd(model.rotation) -sind(model.rotation); sind(model.rotation) cosd(model.rotation)];
        lambda_squared = [(model.lambda(1).^2);(model.lambda(2).^2)];
        inside = sum((rotation_matrix*(pos_1'-pos_2')).^2./lambda_squared).^0.5;
        covariance = exp(-1.*sqrt(5).*inside).*model.variance.*(1+sqrt(5).*inside+5.*inside.^2./(3));
    else
        error('Invalid argument of model.name')
    end
end

function [idx_1,idx_2] = get_box_algorithm(flags,model)
    %Computing the box expicitly is slower due to Matrix Matrix products.
    center = rand(1,2);
    idx_1 = [];
    idx_2 = [];
    for k=1:flags.N_p
        if(abs((flags.pos_all(k,1)/model.size(1))-center(1))<=flags.kappa && abs((flags.pos_all(k,2)/model.size(2))-center(2))<=flags.kappa)       
            idx_1 = [idx_1, k];
        else
            idx_2 = [idx_2, k];
        end
    end
end

function [mu_tilde,covariance_tilde] = compute_mean_and_covariance(idx_1,idx_2,flags)
    % In contrast to the paper, mu_tilde, covariance_tilde do not have an index.
    % Implicitly, both have the subscript value 1.
    if isempty(idx_2)  %No points given, all points are resampled. Hence, the known covarinace can be used.
        covariance_tilde = flags.covariance_matrix;
        mu_tilde = flags.mu;
    else
        % compute mean z and convariance k of not sampled points.
        
        % That implements Kriging (also named gaussian process regression)
        if numel(idx_2) < numel(idx_1) %Decide which calculation method is faster.                    
            temp = flags.covariance_matrix(idx_1,idx_2)/flags.covariance_matrix(idx_2,idx_2); % That is slow
            mu_tilde = flags.mu(idx_1) + temp*(flags.theta(idx_2)-flags.mu(idx_2));
            covariance_tilde = flags.covariance_matrix(idx_1,idx_1)- temp*flags.covariance_matrix(idx_2,idx_1);
        else
            covariance_tilde = inv(flags.covariance_matrix_inverse(idx_1,idx_1)); %linear algebra trick: see Wolfgang Nowak dissertation, Appendix A, formula (A.6)
            mu_tilde = flags.mu(idx_1) - covariance_tilde * flags.covariance_matrix_inverse(idx_1,idx_2)*(flags.theta(idx_2)-flags.mu(idx_2)); %linear algebra trick: see Wolfgang Nowak dissertation, Appendix A, formula (A.7)   
        end

    end
    %Numerical errors need to be filtered out:
    covariance_tilde = 0.5*(covariance_tilde+covariance_tilde'); %Make symetric
    covariance_tilde = max (covariance_tilde,0); %Cut negative vaules which are -1*10^-16
end

function [objective_function] = evaluate_hyper_parameters(flags,model)
    recent_thetas =  flags.theta_recently;
     counter_NaN = 0;
    for i=1:size(flags.theta_all,1)
        %% Estimate efficiency (see: Geyer, C. J. 19992. Practical Markov chain Monte Carlo (with discussion). Statistical Science, 7:473-511)
        acf  = autocorr(recent_thetas(i,:),'NumLags',size(recent_thetas,2)-1);        
        if(mod(size(acf,1),2)==0) %even
                inital_positive_sequence(:) = acf(:,1:2:end-1)+acf(:,2:2:end);
        else %odd
                inital_positive_sequence(:) = acf(1:2:end-2)+acf(2:2:end-1);
        end
        index = find(inital_positive_sequence(:)<0);
        if(~isempty(index))
            inital_positive_sequence(index:end) = 0;
        else
            index =  size(inital_positive_sequence,1);
        end
        x = [1:index(1)]';
        if (size(inital_positive_sequence(1:index(1)),1)==1) %MATLAB 2018a and 2020a have different behaviour.
            P = [x,inital_positive_sequence(1:index(1))'];
        else
            P = [x,inital_positive_sequence(1:index(1))];
        end
        try
            [k,av] = convhull(P);
            max_index_k = find(max(k)==k);
            inital_convex_sequence = interp1(k(1:max_index_k),inital_positive_sequence(k(1:max_index_k)),x);
            efficiency_pre =  sum(inital_convex_sequence);
            efficiency_pre_all(i) = efficiency_pre;
        catch
            efficiency_pre = 0;
            efficiency_pre_all(i) = 0;  %That is the variance = 0 case.
        end

        %% compute mean of Efficiencies
        variance = var(recent_thetas(i,:));
        variance_all(i) = variance;
        if (numel(find(recent_thetas(i,:)==recent_thetas(i,1)))==numel(recent_thetas(i,:))) %All values are the same. Due to numerical errors that does not follow that var = 0;
            counter_NaN = counter_NaN + 1;
            objective_pre(i) = NaN;
        else
            objective_pre(i) = efficiency_pre/ sqrt(variance); %Try efficiency
        end

    end
    if(counter_NaN/size(flags.theta_all,1)>0.0)
        disp(['CARE! id' num2str(model.id) ' A lot of values did not change during the last optimizer step. This can lead to wrong results. Increase "flags.parameter_tuning_frequency" to remedy this error. counter_NaN = ' num2str(counter_NaN)]);               
    end

    %Add worst possible value instead of NaN. 
    objective_pre(isnan(objective_pre))=max(objective_pre);
    objective_function =   1 / (1 + 2 * nanmean(objective_pre));
    if (isnan(objective_function))
        objective_function = 0;
    end
end

%% Check inputs
function [model,solver,flags] = check_inputs(number_of_samples,model,solver,flags)
    if(~isfield(flags,'debug'))
        flags.debug = false;
    end
    if(~isfield(flags,'x_pos'))
        flags.x_pos =  model.size(1)/model.discretization(1)/2:model.size(1)/model.discretization(1):model.size(1);
    end
    if(~isfield(flags,'y_pos'))
        flags.y_pos =  model.size(2)/model.discretization(2)/2:model.size(2)/model.discretization(2):model.size(2);
    end
    model.x_pos = flags.x_pos;
    model.y_pos = flags.y_pos;
    if(~isfield(flags,'pos_all'))
        x_pos =  model.size(1)/model.discretization(1)/2:model.size(1)/model.discretization(1):model.size(1);
        y_pos =  model.size(2)/model.discretization(2)/2:model.size(2)/model.discretization(2):model.size(2);
        [X,Y] = meshgrid(x_pos,y_pos);
        flags.pos_all = [reshape(X,[numel(X),1]),reshape(Y,[numel(Y),1])];
    end
    if(~isfield(flags,'N_p'))
        flags.N_p = size(flags.pos_all,1);
    end
    % Compute covariance matrix
    if(~isfield(flags,'covariance_matrix'))
        flags.covariance_matrix = zeros(flags.N_p);
        for i=1:flags.N_p
            flags.covariance_matrix(i,:) = compute_covariance_matrix(flags.pos_all(i,:),flags.pos_all(:,:),model);
        end 
    end
     if(~isfield(flags,'current_sample_number'))
        flags.current_sample_number=0;
    end
    if(flags.current_sample_number == 0 && isfield(solver,'set_up_solver_function_handle'))
        [model,flags,solver] = solver.set_up_solver_function_handle(model,flags,solver);
        randomize(model)
    end
    
    if(~isfield(flags,'covariance_matrix_inverse'))
        flags.covariance_matrix_inverse = inv(flags.covariance_matrix); % That is slow
    end
    if(~isfield(flags,'mu'))
        flags.mu = ones(flags.N_p,1)*model.mu; % Assume constant mean. That can be changed.
    end
    if(~isfield(flags,'theta'))
        flags.theta = mvnrnd(flags.mu,flags.covariance_matrix)';  %Start with random theta
    end
    if(~isfield(flags,'likelihood'))
        flags.log_likelihood = solver.log_likelihood_function_handle(flags.theta,model,solver,flags); 
    end
    if(~isfield(flags,'theta_all')&&~isfield(flags,'number_of_samples_max'))
        flags.number_of_samples_max = number_of_samples;
    end
    if(~isfield(flags,'theta_all'))
        flags.theta_all=zeros(flags.N_p,floor(flags.number_of_samples_max/model.saving_distance));
    end
    if(~isfield(flags,'log_likelihood_all'))
        flags.log_likelihood_all=zeros(1,floor(flags.number_of_samples_max/model.saving_distance));
    end
    if(~isfield(flags,'number_accepted'))
        flags.number_accepted=0;
    end
    if(~isfield(flags,'parameter_tuning_frequency'))
        flags.parameter_tuning_frequency = 1000;
    end
    if(~isfield(flags,'beta')||~isfield(flags,'kappa'))
        if(isfield(flags,'beta'))
            flags.tested_beta  = flags.beta;
            flags.beta_is_constant = 1;
        else
            flags.tested_beta  = 10^(rand(1)*-2);
            flags.beta_is_constant = 0;
        end
        if(isfield(flags,'kappa'))
            flags.tested_kappa  = flags.kappa;
            flags.kappa_is_constant = 1;
        else
            flags.tested_kappa = 10^(rand(1)*-2);
            flags.kappa_is_constant = 0;
        end
        
        flags.enable_parameter_tuning = true;  
        flags.delta_beta_proportion = 1.414; 
        flags.delta_kappa_proportion = 1.414; 
        flags.max_tuning_iterations = 50; 
        
        %Start with parameter tuning
        flags.kappa = flags.tested_kappa;
        flags.beta  = min(max(flags.tested_beta  * flags.delta_beta_proportion ,0),1); 
        flags.theta_recently = zeros(flags.N_p,flags.parameter_tuning_frequency); 
        
    else
        flags.enable_parameter_tuning = false;  
    end
    if(~isfield(flags,'finished'))
        flags.finished = false;
    end
end

function randomize(model)
    if(isfield(model, 'id'))
        rng(sum(100*clock)*model.id,'Twister')
    else
        rng(sum(100*clock),'Twister')
    end
end