function [model,flags,solver] = set_up_solver_example(model,flags,solver)
%This function gets called once before the start of the algorithm. You can
%setup your simulation tool here. 

if(isfield(solver,'kriging')&&solver.kriging==1)
    %% Kriging
    meas_loc_sorted=[10,10;10,30;10,50;10,70;10,90;30,10;30,30;30,50;30,70;30,90;50,10;50,30;50,50;50,70;50,90;70,10;70,30;70,50;70,70;70,90;...
        90,10;90,30;90,50;90,70;90,90;20,20;20,40;20,60;20,80;40,20;40,40;40,60;40,80;60,20;60,40;60,60;60,80;80,20;80,40;80,60;80,80]*50.5;
    meas_loc = meas_loc_sorted;
    for i=1:size(meas_loc,1)
        id = find(sum((meas_loc(i,:) - flags.pos_all).^2,2)==min(sum((meas_loc(i,:) - flags.pos_all).^2,2))); %Which id is clostest to the measurement point.
        solver.ids_of_measurements(i) = id(1); %Take the first one.  
    end
    
    %Create reference solution
    rng(0)
    load('reference_kriging');
    theta_correct_size = reshape(reference_kriging,[50,50]);
    solver.std         = 0.05;
    solver.measurements = theta_correct_size(solver.ids_of_measurements)'+randn(numel(solver.ids_of_measurements),1).*solver.std;   %add normal noise to measurements
    rng shuffle
    visualize(model,flags,solver,theta_correct_size)
else
    %% FEM-Solver
    % init_FORWARD
    n           = 0.35;
    al          = 2.5;
    at          = 0.25;
    Dm          = 1e-9;

    % init_SOLVERS
    FEMresh     = 1e-10;

    % init_GENERATE_WELLS
    well_loc_x      = [10,70,40,40]*50; %In meter
    well_loc_y      = [47,47,71,21]*50; %In meter
    pump_strength   = [120,70,90,90]; %m^3/d
    pump_strength   = pump_strength/86400*1000; %liter/s

    %% initializing FEM Solver
    % This for parallel
    FEM_path = ['temp' filesep 'FEM_Solver_id' num2str(model.id)];
    %mkdir(FEM_path)
    copyfile('FEM_Solver',FEM_path)

    solver.n = n;
    solver.al = al;
    solver.at = at;
    solver.Dm = Dm;
    solver.FEMresh = FEMresh;
    solver.well_loc_x = well_loc_x;
    solver.well_loc_y = well_loc_y;
    solver.pump_strength = pump_strength;
    solver.FEM_path = FEM_path;

    return_path = pwd;
    cd (FEM_path)
    Fct_define_grid_and_solver(solver, model);
    cd (return_path)

    [model,flags,solver] =  get_reference_solution(model,flags,solver);
end
end


function [model,flags,solver] = get_reference_solution(model,flags,solver)
    solver.benchmark_type = 'mild';
    if(strcmp(solver.benchmark_type,'matern'))
        field = load(['FEM_Solver' filesep 'reference_matern.mat']);
        theta_correct_size = field.reference_matern;
    else
        load(['FEM_Solver' filesep 'model_reference_large.lpf']);
        load(['FEM_Solver' filesep 'model_reference_mild.lpf']);
        if(strcmp(solver.benchmark_type,'mild'))
            theta_temp = flipud(log(model_reference_mild));
        elseif(strcmp(solver.benchmark_type,'large'))
            theta_temp = flipud(log(model_reference_large));
        else
            error('benchmark_type can either be mild or large. Please choose which test case you want to use!')
        end

        if(model.discretization(1)== 50 && model.discretization(2)== 50)
            theta_temp = interp2(theta_temp); 
            theta_correct_size = theta_temp(2:4:end,2:4:end); %Only works for 50 by 50.    
        elseif(model.discretization(1)== 100 && model.discretization(2)== 100)
            theta_correct_size = theta_temp;

        else
            error('This example only works for model.discretization = [50,50] or model.discretization = [100,100].');
        end

    end
    
    
    %% calculate refrence hydraulic heads with FEM

    %call FEM simulation
    return_path = pwd;
    cd (solver.FEM_path)
    [h_clean,drawdown] = Fct_FEM_interface(theta_correct_size,solver);
    cd(return_path);

    %calculate head in cell centers based on simulation results
    h = h_clean + drawdown;

    h = reshape(h,[model.discretization(2)+1,model.discretization(1)+1]);
    h = interp2(h); 
    h = h(2:2:end,2:2:end);
    h = reshape(h,[prod(model.discretization),1]);
    
    %% Save head at measurement locations
   
    % Define Measurements
    solver.std = 0.05;
    model.meas = 41;
    meas_loc_random=[70,29;23,80;45,14;12,86;28,16;67,46;69,78;23,89;63,22;12,28;64,63;38,60;6,13;61,38;10,24;8,56;20,71;92,82;37,34;71,36; ...
        88,59;9,34;74,88;46,47;82,29;8,92;24,54;58,65;34,72;84,12;5,94;32,88;50,37;7,30;17,20;40,28;21,7;50,95;68,93;49,16;69,57]*50.5; %In meters
    meas_loc_sorted=[10,10;10,30;10,50;10,70;10,90;30,10;30,30;30,50;30,70;30,90;50,10;50,30;50,50;50,70;50,90;70,10;70,30;70,50;70,70;70,90;...
        90,10;90,30;90,50;90,70;90,90;20,20;20,40;20,60;20,80;40,20;40,40;40,60;40,80;60,20;60,40;60,60;60,80;80,20;80,40;80,60;80,80]*50.5;
 
    if(model.meas==41)
        meas_loc = meas_loc_sorted;
    elseif(model.meas==16)
        meas_loc = [20,20;20,40;20,60;20,80;40,20;40,40;40,60;40,80;60,20;60,40;60,60;60,80;80,20;80,40;80,60;80,80]*50.5;
    elseif(model.meas==4)
        meas_loc = [30,30;30,70;70,30;70,70]*50.5;
    elseif(model.meas==1)
        meas_loc = [50,50]*50.5;
    elseif(model.meas==61)
        meas_loc = [meas_loc_sorted ;[0,0;0,20;0,40;0,60;0,80;0,100;  100,0;100,20;100,40;100,60;100,80;100,100;   20,0;40,0;60,0;80,0; 20,100;40,100;60,100;80,100 ]*50.5];
    end
    
    
    for i=1:size(meas_loc,1)
        id = find(sum((meas_loc(i,:) - flags.pos_all).^2,2)==min(sum((meas_loc(i,:) - flags.pos_all).^2,2))); %Which id is clostest to the measurement point.
        solver.ids_of_measurements(i) = id(1); %Take the first one.  
    end
    
    rng(0) %Always add the same noise
    if(isfield(solver,'kriging')&&solver.kriging==1)
         solver.measurements = theta_correct_size(solver.ids_of_measurements)'+randn(numel(solver.ids_of_measurements),1).*solver.std;   %add normal noise to measurements
    else
        solver.measurements = h(solver.ids_of_measurements)+randn(numel(solver.ids_of_measurements),1).*solver.std;   %add normal noise to measurements   
    end
    rng shuffle    
    visualize(model,flags,solver,theta_correct_size)
end

function visualize(model,flags,solver,theta_correct_size)
    if(model.visualization_distance>0)
        figure(1)
        subplot(1,2,1)        
        imagesc(model.x_pos,model.y_pos,reshape(theta_correct_size,[model.discretization(2),model.discretization(1)]));set(gca,'YDir','normal');daspect([1 1 1]);colorbar;
        xlabel('x [m]'); ylabel('y [m]'); title('theta');
        subplot(1,2,2); solver.log_likelihood_function_handle(theta_correct_size,model,solver,flags,1); %Inefficient because the function gets called again here.
       try
            suptitle('Reference solution');
        catch
            try
                sgtitle('Reference solution');
            catch
            end
        end
        drawnow;
    end
    
end
