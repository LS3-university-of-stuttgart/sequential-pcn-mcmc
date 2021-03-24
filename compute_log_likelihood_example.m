function [log_likelihood] = compute_log_likelihood_example(theta,model,solver,flags,enable_plot)
%This function gets called every time the likelihood is evaluated. 
%Try to keep this funciton as efficient as possible because it gets
%called thousands of times. 
if(~exist('enable_plot'))
    enable_plot = 0;
end


if(isfield(solver,'kriging')&&solver.kriging==1)
    log_likelihood = sum( log(   1./(sqrt(2*pi*solver.std.^2)))   -         (reshape(theta(solver.ids_of_measurements),[numel(theta(solver.ids_of_measurements)),1])-solver.measurements).^2/(2.*solver.std.^2));
    h = theta;
else

    %% calculate hydraulic heads with FEM

    %call FEM simulation
    return_path = pwd;
    cd (solver.FEM_path)
    [h_clean,drawdown] = Fct_FEM_interface(theta,solver);
    cd(return_path);

    %calculate head in cell centers based on simulation results
    h = h_clean + drawdown;

    h = reshape(h,[model.discretization(2)+1,model.discretization(1)+1]);
    h2 = imresize(h,[model.discretization(2),model.discretization(1)],'bilinear');
    h = interp2(h); 
    h = h(2:2:end,2:2:end);
    h = reshape(h,[prod(model.discretization),1]);

    %calculate log likelihood.
    log_likelihood = sum( log(   1./(sqrt(2*pi*solver.std.^2)))   -         (h(solver.ids_of_measurements)-solver.measurements).^2/(2.*solver.std.^2));
end

%plot
if enable_plot
    imagesc(model.x_pos,model.y_pos,reshape(h,[model.discretization(2),model.discretization(1)]));
    set(gca,'YDir','normal');daspect([1 1 1]);colorbar;
    xlabel('x [m]'); ylabel('y [m]');
    if(isfield(solver,'kriging')&&solver.kriging==1)
        title('Hydraulic conductivity');
    else
        title('Hydraulic head');
    end
    hold on
    plot(flags.pos_all(solver.ids_of_measurements,1),flags.pos_all(solver.ids_of_measurements,2),'ok','LineWidth',2,...
                                        'MarkerSize',5,...
                                        'MarkerEdgeColor','k',...
                                        'MarkerFaceColor','w');
    hold off
end

end