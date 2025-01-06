%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_plot_solution.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIESIM_plot_solution_2D(solution, psize, uinput, grid, opts);
% This routine outputs and plots solution of ODE and PDE states in 2D


% Inputs: 
% 1) solution 
% solution is a structure with the following fields
% --- solution.tf - scalar - actual final time of the solution
% --- solution.final.pde{1,2} - array containing all pde solutions at a final time
% --- solution.final.pde{1} - array containing the solution for states that are only the functions of one variable - 
% it is array of size (N+1) x (nx+ny), nx - number of states depending only on s1, 
%                                      ny - number of states depending only on s2 
% --- solution.final.pde{2} - array containing the solution for states that are the functions of two variables - 
% it is array of size (N+1) x (N+1) x ns, ns - number of states depending on both s1 and s2
% --- solution.final.ode - array of size no - ode solution at a final time 
% --- solution.final.observed - array of size noo  - final value of observed outputs
% --- solution.final.regulated - array of size nro  - final value of regulated outputs

% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.pde{1,2} - array containing all time-dependent pde solutions 
% --- solution.timedep.pde{1} - array containing the solution for states that are only the functions of one variable - 
% it is array of size (N+1) x (nx+ny) x Nsteps, nx - number of states depending only on s1, 
%                                      ny - number of states depending only on s2 
% --- solution.timedep.pde{2} - array containing the solution for states that are the functions of two variables - 
% it is array of size (N+1) x (N+1) x ns x Nsteps, ns - number of states depending on both s1 and s2
% --- solution.timedep.ode - array of size no x Nsteps - time-dependent solution of no ODE states
% --- solution.timedep.observed - array of size noo x Nsteps -
%     time-dependent value of observed outputs
% --- solution.timedep.regulated - array of size nro x Nsteps -
%     time-dependent value of regulated outputs

% 2) psize - all variables defining the size of the PIE problem
% 3) uinput -  user-defined boundary inputs, forcing functions and initial conditions
% 4) grid - physical and computational grid for states differentiable up to order zero (corresponding to a orimary = PDE state discretization)
% 5) opts - options for simulation parameters
% 
% Outputs: none


% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% DJ, 07/17/2024 - Bugfix for extraction of exact solutions.
% DJ, 01/01/2025 - Change plot specifications;

syms sx sy;

line_width = 2.5;
marker_size = 2*line_width;

   % Output solution of ODE states

for i=1:psize.no
formatSpec = 'Solution of an ODE state %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.ode(i));
end

% Output values for regulated outputs
for i=1:psize.nro
formatSpec = 'Value of regulated output %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.regulated(i));
end

% Output values for observed outputs
for i=1:psize.noo
formatSpec = 'Value of observed output %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.observed(i));
end

% Don't plot if not requested
if ~strcmp(opts.plot,'yes')
    return
end
% Temporal evolution of solution is only available with BDF scheme
if opts.intScheme~=1
    disp("Temporal evolution of solution is only available when using the BDF scheme; only final states are plotted.")
elseif opts.tf==0
    disp("Final time is t=0; no plots produced.")
    return
else
    dtime=solution.timedep.dtime;
    Nplot = min(length(dtime),100);  % plot at at most 100 time points.
    t_idcs = floor(linspace(1,length(dtime),Nplot));
    tval = dtime(t_idcs);
end
figs = {};

n_pde_tot = sum([sum(psize.nx),sum(psize.ny),sum(psize.n)]);
if opts.intScheme==1 && psize.no>0 && ~isempty(solution.timedep.ode)
    % Plot evolution of ODE states in a single plot
    odesol = solution.timedep.ode';
    figure('Position',[200 150 800 400]);
    set(gcf, 'Color', 'w');
    box on
    hold on
    for ii=1:psize.no
        X_ii = odesol(t_idcs,ii)';
        plot(tval,X_ii,'DisplayName',['$x_',num2str(ii),'(t)$'],'LineWidth',line_width);
    end
    hold off
    %plot(dtime,odesol','-S','LineWidth',2,'DisplayName',['$x_',num2str(i),'$(t)']); 
    if strcmp(opts.type,'DDE')
        title('DDE State Evolution','FontSize',16,'Interpreter','latex');
    else
        title('ODE State Evolution','FontSize',16,'Interpreter','latex');
    end
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$x$','FontSize',15,'Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    % ax = gca;   ax.FontSize = 16;
    % H=gca;
    % H.LineWidth=3;
    if psize.no>1
        [leg,hobj] = legend('FontSize',15,'Interpreter','latex');
        %set(leg,'Box','off','Location','north','FontSize',16);
        %hl = findobj(hobj,'type','line');
        %set(hl,'LineWidth',2);  
    end
elseif opts.intScheme==1 && psize.no>0
    disp('ODE solution is infinite. Unable to plot.')
end

% % Plot the observed outputs evolution in a single plot.
colors_y = [0, 0, 1;
            0, 0.5, 0;
            1, 0, 0;
            0, 0.75, 0.75;
            0.75, 0, 0.75;
            0.75, 0.75, 0;
            0.25, 0.25, 0.25];
if opts.intScheme==1 && (psize.noo>0) && ~isempty(solution.timedep.observed)
    fig2 = figure('Position',[200 150 800 400]);
    box on
    title('Observed Output Evolution','FontSize',16,'Interpreter','latex');
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$y$','FontSize',15,'Interpreter','latex');
    hold on
    for ii=1:psize.noo
        y_ii = solution.timedep.observed(ii,t_idcs);
        if ii<=size(colors_y,1)
            plot(tval,y_ii,'Color',colors_y(psize.noo-ii+1,:),'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',line_width);   
        else
            plot(tval,y_ii,'DisplayName',['$y_',num2str(ii),'(t)$'],'LineWidth',line_width);
        end
    end
    hold off
    if psize.noo>1
        legend('FontSize',15,'Interpreter','latex');
    end
    set(gca,'XLim',[min(tval),max(tval)]);
    set(gca,'TickLabelInterpreter','latex');
    set(gcf, 'Color', 'w');
    figs = [figs,{fig2}];

    % y = solution.timedep.observed';
    % for i=1:psize.noo
    %     labels{i} = ['y_',num2str(i)];
    % end            
    % figure;
    % plot(dtime,y,'-S','linewidth',2); 
    % title('Time evolution of observed outputs');
    % xlabel('Time');
    % ylabel('Value of an observed output');
    % ax = gca;
    % ax.FontSize = 16;
    % H=gca;
    % H.LineWidth=3;
    % [leg, hobj]=legend(labels);
    % set(leg,'Box','off','Location','north','FontSize',16);
    % hl = findobj(hobj,'type','line');
    % set(hl,'LineWidth',2);  
    % clear labels;
elseif opts.intScheme==1 && (psize.noo>0)
    disp('Observed output is infinite. Unable to plot.')
end

% % Plot the regulated outputs evolution in a single plot.
colors_z = {'b','g','m','r','k','c','r','y'};
if opts.intScheme==1 && (psize.nro>0) && ~isempty(solution.timedep.regulated)
    fig3 = figure('Position',[200 150 800 400]);
    box on
    title('Regulated Output Evolution','FontSize',16,'Interpreter','latex');
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$z$','FontSize',15,'Interpreter','latex');
    hold on
    for ii=1:psize.nro
        z_ii = solution.timedep.regulated(ii,t_idcs);
        if ii<=length(colors_z)
            plot(tval,z_ii,[colors_z{ii},'-'],'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',line_width);   
        else
            plot(tval,z_ii,'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',line_width);
        end
    end
    hold off
    if psize.nro>1
        legend('FontSize',15,'Interpreter','latex');
    end
    set(gca,'XLim',[min(tval),max(tval)]);
    set(gca,'TickLabelInterpreter','latex');
    set(gcf, 'Color', 'w');
    figs = [figs,{fig3}];

    % z=solution.timedep.regulated';
    % for i=1:psize.nro
    %     labels{i} = ['z_',num2str(i)];
    % end 
    %  figure;
    %  plot(dtime,z,'-S','linewidth',2); 
    %  title('Time evolution of regulated outputs');
    %  xlabel('Time');
    %  ylabel('Value of a regulated output');
    %  ax = gca;
    %  ax.FontSize = 16;
    %  H=gca;
    %  H.LineWidth=3;
    %  [leg, hobj]=legend(labels);
    %  set(leg,'Box','off','Location','north','FontSize',16);
    %  hl = findobj(hobj,'type','line');
    %  set(hl,'LineWidth',2);  
    %  clear labels;
elseif opts.intScheme==1 && (psize.nro>0)
    disp('Regulated output is infinite. Unable to plot.')
end


% % Plot final PDE states, if possible
if strcmp(opts.type,'DDE')
    return
elseif isempty(solution.final.pde) && n_pde_tot>0
    disp("Simulated PDE solution is infinite; no plot produced.")
    return
end
% Declare colors for plots
colors_PDE = [0 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.9290 0.6940 0.1250;
            0.4940 0.1840 0.5560;
            0.4660 0.6740 0.1880;
            0.3010 0.7450 0.9330;
            0.6350 0.0780 0.1840];
if (uinput.ifexact==true)
    Nplot_space=101;
    a = uinput.dom(1,1);    b = uinput.dom(1,2);
    exact_grid_x = linspace(a,b,round(Nplot_space));
    c = uinput.dom(2,1);    d = uinput.dom(2,2);
    exact_grid_y = linspace(c,d,round(Nplot_space));
end

% Keep track of total number of states for indexing
ns_tot = psize.no;

% Plot numerical solution using only markers if exact solution is available
if uinput.ifexact
    line_style = {'d','LineWidth',marker_size};
else
    line_style = {'-d','MarkerSize',marker_size,'LineWidth',line_width};
end

% Plot 1D states, differentiable in x
if sum(psize.nx)>0
    ns = sum(psize.nx);
    fig_width = min(ns*600,1800);
    fig4 = figure('Position',[200 150 fig_width 400]);
    set(gcf, 'Color', 'w');
    box on
    if ns>1
        sgtitle('Simulated Final 1D PDE State','Interpreter','latex','FontSize',16)
    end
    for n=1:ns
        subplot(1,ns,n);
        box on
        if n+ns_tot<=size(colors_PDE,1)
            plot(grid.phys(:,1),solution.final.pde{1}(:,n),line_style{:},'Color',colors_PDE(n+ns_tot,:),'DisplayName','Numerical solution');
        else
            plot(grid.phys(:,1),solution.final.pde{1}(:,n),line_style{:},'DisplayName','Numerical solution');
        end
        if uinput.ifexact
            hold on
            exsol_grid_time = subs(uinput.exact(n+ns_tot),sx,exact_grid_x);
            exsol_grid_x = double(subs(exsol_grid_time,solution.tf));
            plot(exact_grid_x,exsol_grid_x,'k-','LineWidth',line_width,'DisplayName','Analytic solution'); 
            legend('FontSize',13,'Interpreter','latex');
            hold off
        end
        xlabel('$s_{1}$','FontSize',15,'Interpreter','latex');
        ylabel('$\mathbf{x}$','FontSize',15,'Interpreter','latex');
        if ns==1 && n_pde_tot==1
            title(['Simulated Final 1D PDE State, $\mathbf{x}(t=',num2str(solution.timedep.dtime(end)),',s_{1})$'],'FontSize',15,'Interpreter','latex');
        elseif ns==1
            title(['Simulated Final 1D PDE State, $\mathbf{x}_',num2str(n+ns_tot),'(t=',num2str(solution.timedep.dtime(end)),',s_{1})$'],'FontSize',15,'Interpreter','latex');
        else
            title(['$\mathbf{x}_',num2str(n+ns_tot),'(t=',num2str(solution.timedep.dtime(end)),',s_{1})$'],'FontSize',15,'Interpreter','latex');
        end
        set(gca,'XLim',[min(grid.phys(:,1)),max(grid.phys(:,1))]);
        set(gca,'TickLabelInterpreter','latex');
    end
    figs = [figs,{fig4}];
    ns_tot = ns_tot+ns;
end

% Plot 1D states, differentiable in y
if sum(psize.ny)>0
    ns = sum(psize.ny);
    fig_width = min(ns*600,1800);
    fig5 = figure('Position',[200 150 fig_width 400]);
    set(gcf, 'Color', 'w');
    box on
    if ns>1
        sgtitle('Simulated Final 1D PDE State','Interpreter','latex','FontSize',16)
    end
    for n=1:ns
        subplot(1,ns,n);
        box on
        if sum(psize.nx)+n<=size(colors_PDE,1)
            plot(grid.phys(:,2),solution.final.pde{1}(:,sum(psize.nx)+n),line_style{:},'Color',colors_PDE(n+ns_tot,:),'DisplayName','Numerical solution');
        else
            plot(grid.phys(:,2),solution.final.pde{1}(:,sum(psize.nx)+n),line_style{:},'DisplayName','Numerical solution');
        end
        if uinput.ifexact
            hold on
            exsol_grid_time = subs(uinput.exact(n+ns_tot),sy,exact_grid_y);
            exsol_grid_y = double(subs(exsol_grid_time,solution.tf));
            plot(exact_grid_y,exsol_grid_y,'k-','LineWidth',line_width,'DisplayName','Analytic solution'); 
            legend('FontSize',13,'Interpreter','latex');
            hold off
        end
        xlabel('$s_{2}$','FontSize',15,'Interpreter','latex');
        ylabel('$\mathbf{x}$','FontSize',15,'Interpreter','latex');
        if ns==1 && n_pde_tot==1
            title(['Simulated Final 1D PDE State, $\mathbf{x}(t=',num2str(solution.timedep.dtime(end)),',s_{2})$'],'FontSize',15,'Interpreter','latex');
        elseif ns==1
            title(['Simulated Final 1D PDE State, $\mathbf{x}_',num2str(n+ns_tot),'(t=',num2str(solution.timedep.dtime(end)),',s_{2})$'],'FontSize',15,'Interpreter','latex');
        else
            title(['$\mathbf{x}_',num2str(n+ns_tot),'(t=',num2str(solution.timedep.dtime(end)),',s_{2})$'],'FontSize',15,'Interpreter','latex');
        end
        set(gca,'XLim',[min(grid.phys(:,2)),max(grid.phys(:,2))]);
        set(gca,'TickLabelInterpreter','latex');
    end
    figs = [figs,{fig5}];
    ns_tot = ns_tot+ns;
end

% Plot final 2D states using isosurface
if sum(psize.n)>0
    ns = sum(psize.n);

    % PIESIM stores nD solutions by columns, i.e., e.g., for u(x,y), the
    % entries with fixed y and varying x are stored in columns, and the
    % entries with fixed x and varying y are stored in rows

    % MATLAb surf routine does the opposite, so we plot the transpose of
    % the solution vector to get a correct orienation in space 

    % Plot the numerical value
    fig_width = min(ns*600,1800);
    fig6 = figure('Position',[200 150 fig_width 400]);
    if ns>1
        sgtitle('Simulated Final 2D PDE State','Interpreter','latex','FontSize',16)
    end
    set(gcf, 'Color', 'w');
    box on
    for n=1:ns
        ax1 = subplot(1,ns,n,'Parent',fig6);
        surf(ax1,grid.phys(:,1),grid.phys(:,2),solution.final.pde{2}(:,:,n)','FaceAlpha',0.75,'Linestyle','--','FaceColor','interp');
        h = colorbar(ax1);
        %colormap jet
        xlabel(ax1,'$s_{1}$','FontSize',15,'Interpreter','latex');
        ylabel(ax1,'$s_{2}$','FontSize',15,'Interpreter','latex');
        zlabel(ax1,'$\mathbf{x}$','FontSize',15,'Interpreter','latex');
        box on
        if ns==1 && n_pde_tot==1
            title(ax1,['Simulated Final PDE State $\mathbf{x}(',num2str(opts.tf),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
        elseif ns==1
            title(ax1,['Simulated Final 2D PDE State $\mathbf{x}_',num2str(n)+ns_tot,'(',num2str(opts.tf),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
        else
            title(ax1,['$\mathbf{x}_',num2str(n)+ns_tot,'(',num2str(opts.tf),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
        end
        ax1.XLim = [min(grid.phys(:,1)),max(grid.phys(:,1))];
        ax1.YLim = [min(grid.phys(:,2)),max(grid.phys(:,2))];
        ax1.TickLabelInterpreter = 'latex';
    end
    figs = [figs,{fig6}];

    % Plot the analytic PDE solution
    if uinput.ifexact
        fig7 = figure('Position',[200 150 fig_width 400]);
        if ns>1
            sgtitle('Exact Final 2D PDE State','Interpreter','latex','FontSize',16)
        end
        set(gcf, 'Color', 'w');
        box on
        fig8 = figure('Position',[200 150 fig_width 400]);
        if ns>1
            sgtitle('Error in Final 2D PDE State','Interpreter','latex','FontSize',16)
        end
        set(gcf, 'Color', 'w');
        box on
        for n=1:ns
            % Compute the value of the exact solution at final time
            exsol_grid_time = subs(subs(uinput.exact(n+ns_tot),sx,exact_grid_x'),sy,exact_grid_y);
            exsol_grid = double(subs(exsol_grid_time,solution.tf));
            exsol_numgrid_time = subs(subs(uinput.exact(n+ns_tot),sx,grid.phys(:,1)),sy,grid.phys(:,2)');
            exsol_numgrid = double(subs(exsol_numgrid_time,solution.tf));
            
            % Plot the exact solution
            ax2 = subplot(1,ns,n,'Parent',fig7);
            surf(ax2,exact_grid_x,exact_grid_y,exsol_grid','FaceAlpha',0.75,'FaceColor','interp');
            h = colorbar(ax2);
            %colormap jet
            xlabel(ax2,'$s_{1}$','FontSize',15,'Interpreter','latex');
            ylabel(ax2,'$s_{2}$','FontSize',15,'Interpreter','latex');
            zlabel(ax2,'$\mathbf{x}_{true}$','FontSize',15,'Interpreter','latex');
            box on
            if ns==1 && n_pde_tot==1
                title(ax2,['Exact Final PDE State $\mathbf{x}_{true}(',num2str(opts.tf),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
            elseif ns==1
                title(ax2,['Exact Final 2D PDE State $\mathbf{x}_{true,',num2str(n)+ns_tot,'}(',num2str(opts.tf),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
            else
                title(ax2,['$\mathbf{x}_{true,',num2str(n)+ns_tot,'}(',num2str(opts.tf),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
            end
            ax2.XLim = [min(grid.phys(:,1)),max(grid.phys(:,1))];
            ax2.YLim = [min(grid.phys(:,2)),max(grid.phys(:,2))];
            ax2.TickLabelInterpreter = 'latex';

            % Also plot the error in the numerical solution
            ax3 = subplot(1,ns,n,'Parent',fig8);
            surf(ax3,grid.phys(:,1),grid.phys(:,2),exsol_numgrid'-solution.final.pde{2}(:,:,n)','FaceAlpha',0.75,'Linestyle','--','FaceColor','interp');
            h = colorbar(ax3);
            %colormap jet
            xlabel(ax3,'$s_{1}$','FontSize',15,'Interpreter','latex');
            ylabel(ax3,'$s_{2}$','FontSize',15,'Interpreter','latex');
            zlabel(ax3,'$\mathbf{x}_{true}-\mathbf{x}_{num}$','FontSize',15,'Interpreter','latex');
            box on
            if ns==1 && n_pde_tot==1
                title(ax3,['Error in Final PDE State $\mathbf{x}(',num2str(opts.tf),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
            elseif ns==1
                title(ax3,['Error in Final 2D PDE State $\mathbf{x}_{',num2str(n)+ns_tot,'}(',num2str(opts.tf),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
            else
                title(ax3,['$\mathbf{x}_{true,',num2str(n)+ns_tot,'}(',num2str(opts.tf),',s_{1},s_{2})-\mathbf{x}_{',num2str(n)+ns_tot,'}(',num2str(opts.tf),',s_{1},s_{2})$'],'FontSize',15,'Interpreter','latex');
            end
            ax3.XLim = [min(grid.phys(:,1)),max(grid.phys(:,1))];
            ax3.YLim = [min(grid.phys(:,2)),max(grid.phys(:,2))];
            ax3.TickLabelInterpreter = 'latex';

            figs = [figs,{fig7},{fig8}];

        end            
    end
end

end
