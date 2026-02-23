%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_plot_solution.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figs = PIESIM_plot_solution(solution, psize, uinput, grid, opts);
% This routine outputs and plots solution of ODE and PDE states

% Inputs: 
% 1) solution 
% solution is a strucutre with the following fields
% --- solution.tf - scalar - actual final time of the solution
% --- solution.final.primary{1,2} - cell array containing final value of
% solution states
% --- solution.final.primary{1} - array of size n0 - ode (finite-dimensional) solutions at a final time 
% --- solution.final.primary{2} - array of size (N+1) x nx - pde (distributed state) solutions at a final time
% --- solution.final.observed{1,2} - cell array containing final value of observed outputs 
% --- solution.final.observed[1} - array of size no0 containing final
%      value of finite-dimensional observed outputs
% --- solution.final.observed{2} - array of size (N+1) x nox 
%      containing final value of infinite-dimensional observed outputs
% --- solution.final.regulated{1,2} - cell array containing final value of regulated outputs 
% --- solution.final.regulated[1} - array of size nr0 containing
%     final value of finite-dimensional regulated outputs
% --- solution.final.regulated{2} - array of size (N+1) x nrx 
%      containing final value of infinite-dimensional regulated ooutputs


% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.primary{1,2} - cell array containing final value of
% solution states
% --- solution.timedep.primary{1} - array of size no x Nsteps - time-dependent solutions of no ODE (finite-dimensional) states
% --- solution.timedep.primary{2} - array of size (N+1) x ns x Nsteps - time-dependent
%     solution of ns PDE (distributed) states of the primary PDE system
% --- solution.timedep.observed{1,2} - cell array containing time-dependent observed outputs 
% --- solution.timedep.observed[1} - array of size no0 x Nsteps -
%     time-dependent value of finite-dimensional observed outputs
% --- solution.timedep.observed{2} - array of size (N+1) x nox x Nsteps
%      containing infinite-dimensional observed outputs
% --- solution.timedep.regulated{1,2} - cell array containing time-dependent regulated outputs 
% --- solution.timedep.regulated[1} - array of size nr0 x Nsteps -
%     time-dependent value of finite-dimensional regulated outputs
% --- solution.timedep.regulated{2} - array of size (N+1) x nrx x Nsteps
% containing infinite-dimensional regulated ooutputs

% 2) psize - size of the PIE problem: nw, nu, nf, n0
% 3) uinput -  user-defined boundary inputs, forcing functions and initial conditions
% 4) grid - physical grid for primary states
% 5) opts - options for simulation parameters
% 
% Outputs: none


% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_19_2021
% YP, 6/16/2022 - added outputs for observed and regulated outputs; 
% DJ, 12/26/2024 - Plot PDE state evolution using surf;
% DJ, 01/01/2025 - Change plot specifications;
% YP, 06/02/2025 - Made modifications so that final plots are produced at a
% zero final time (tf=0);
% YP, 01/02/2026 - Updated plotting to distinguish between finite- and
% infinite-dimensinal outputs;
% YP, 02/17/2026 - Renamed solution.final.ode/solution.final.pde into
% solution.final.primary{1,2}; renamed solution.timedep.ode/solution.timedep.pde into
% solution.timedep.primary{1,2}

syms sx;

line_width = 2.5;
marker_size = 2*line_width;

   
% % Display final value of ODE state, if possible.
if(psize.n0>0)
    if (~isempty(solution.final.primary{1}))
        for i=1:psize.n0
            formatSpec = 'Solution of ODE state %s at a final time %f is %8.4f\n';
            fprintf(formatSpec,num2str(i),solution.tf, solution.final.primary{1}(i));
        end
        plot_ode = true;
    else
        disp('ODE solution is infinite. Unable to plot. Numerical value is not returned.')
        plot_ode = false;
    end
end

% % Display final value of regulated output, if possible.
if(psize.nr0>0)
    if (~isempty(solution.final.regulated{1}))
        for i=1:psize.nr0
            formatSpec = 'Value of regulated finite-dimensional output %s at a final time %f is %8.4f\n';
            fprintf(formatSpec,num2str(i),solution.tf, solution.final.regulated{1}(i));
        end
        plot_z = true;
    else
        disp('Regulated output is infinite. Unable to plot. Numerical value is not returned.')
        plot_z = false;
    end
end

% % Display final value of observed output, if possible.
if(psize.no0>0)
    if (~isempty(solution.final.observed{1}))
        for i=1:psize.no0
            formatSpec = 'Value of observed finite-dimensional output %s at a final time %f is %8.4f\n';
            fprintf(formatSpec,num2str(i),solution.tf, solution.final.observed{1}(i));
        end
        plot_y = true;
    else
        disp('Observed output is infinite. Unable to plot. Numerical value is not returned.')
        plot_y = false;
    end
end

% Don't plot if not requested
if ~strcmp(opts.plot,'yes')
    return
end
% Temporal evolution of solution is only available with BDF scheme
if opts.intScheme~=1
    disp("Temporal evolution of solution is only available when using the BDF scheme; only final states are plotted.")
elseif opts.tf==0
    disp("Temporal evolution of solution is only available when final time is not equal to zero; only final states are plotted.")
    plot_ode=false;
    plot_y=false;
    plot_z=false;
else
    % Determine time indices at which to plot solution.
    dtime = solution.timedep.dtime;
    Nplot = min(length(dtime),100);  % plot at at most 100 time points.
    t_idcs = floor(linspace(1,length(dtime),Nplot));
    tval = dtime(t_idcs);
end
figs = {};


% % Plot the ODE states evolution in a single plot
if (psize.n0>0) && plot_ode && opts.intScheme==1
    odesol = solution.timedep.primary{1}';

    fig1 = figure('Position',[200 150 800 450]);
    set(gcf, 'Color', 'w');
    box on
    
    hold on
    for ii=1:psize.n0
        X_ii = odesol(t_idcs,ii)';
        plot(tval,X_ii,'DisplayName',['$x_',num2str(ii),'(t)$'],'LineWidth',line_width);
    end
    hold off
    %plot(dtime,odesol,'-S','linewidth',2); 
    if strcmp(opts.type,'DDE')
        title('DDE State Evolution','FontSize',16,'Interpreter','latex');
    else
        title('ODE State Evolution','FontSize',16,'Interpreter','latex');
    end
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$x$','FontSize',15,'Interpreter','latex');
    if psize.n0>1
        legend('FontSize',15,'Interpreter','latex');
    end
    set(gca,'TickLabelInterpreter','latex');
    %ax = gca;
    %ax.FontSize = 16;
    %H = gca;
    %H.LineWidth = 3;
    %[leg, hobj] = legend();
    %set(leg,'Box','off','Location','north','FontSize',16);
    %hl = findobj(hobj,'type','line');
    %set(hl,'LineWidth',2);

    figs = [figs,{fig1}];
end


% % Plot the observed outputs evolution in a single plot.
colors_y = [0, 0, 1;
            0, 0.5, 0;
            1, 0, 0;
            0, 0.75, 0.75;
            0.75, 0, 0.75;
            0.75, 0.75, 0;
            0.25, 0.25, 0.25];
if (psize.no0>0) && plot_y && opts.intScheme==1

    fig2 = figure('Position',[200 150 800 450]);
    box on
    title('Observed Output Evolution','FontSize',16,'Interpreter','latex');
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$y$','FontSize',15,'Interpreter','latex');
    hold on
    for ii=1:psize.no0
        y_ii = solution.timedep.observed{1}(ii,t_idcs);
        if ii<=size(colors_y,1)
            plot(tval,y_ii,'Color',colors_y(psize.no0-ii+1,:),'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',line_width);   
        else
            plot(tval,y_ii,'DisplayName',['$y_',num2str(ii),'(t)$'],'LineWidth',line_width);
        end
    end
    hold off
    if psize.no0>1
        legend('FontSize',15,'Interpreter','latex');
    end
    set(gca,'XLim',[min(tval),max(tval)]);
    set(gca,'TickLabelInterpreter','latex');
    set(gcf, 'Color', 'w');
    figs = [figs,{fig2}];
end

% % Plot the regulated outputs evolution in a single plot.
colors_z = {'b','g','m','r','k','c','r','y'};
if (psize.nr0>0) && plot_z && opts.intScheme==1

    fig3 = figure('Position',[200 150 800 450]);
    box on
    title('Regulated Output Evolution','FontSize',16,'Interpreter','latex');
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$z$','FontSize',15,'Interpreter','latex');
    hold on
    for ii=1:psize.nr0
        z_ii = solution.timedep.regulated{1}(ii,t_idcs);
        if ii<=length(colors_z)
            plot(tval,z_ii,[colors_z{ii},'-'],'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',line_width);   
        else
            plot(tval,z_ii,'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',line_width);
        end
    end
    hold off
    if psize.nr0>1
        legend('FontSize',15,'Interpreter','latex');
    end
    set(gca,'XLim',[min(tval),max(tval)]);
    set(gca,'TickLabelInterpreter','latex');
    set(gcf, 'Color', 'w');
    figs = [figs,{fig3}];
end

% % Plot the PDE states
if ~strcmp(opts.type,'DDE') && ~isempty(solution.final.primary{2})

    colors_PDE = [0 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.9290 0.6940 0.1250;
            0.4940 0.1840 0.5560;
            0.4660 0.6740 0.1880;
            0.3010 0.7450 0.9330;
            0.6350 0.0780 0.1840];

    % Initialize grid points at which to compute exact solution.
    if (opts.ifexact==true)
        Nplot_space=101;
        a = uinput.a;
        b = uinput.b;
        exact_grid = linspace(a,b,round(Nplot_space));
        exsol_grid = double.empty(round(Nplot_space),0);
    end
    nx = sum(psize.n);

    % Plot the PDE states at the final time
    fig_width = min(nx*600,1800);
    fig4 = figure('Position',[200 150 fig_width 450]);
    set(gcf, 'Color', 'w');
    box on
    if nx>1
        sgtitle('PDE State at Final Time','Interpreter','latex','FontSize',16)
    end
    for n=1:nx
        if (opts.ifexact==true)
            exsol_grid_time = subs(uinput.exact(n),sx,exact_grid);
            exsol_grid = double(subs(exsol_grid_time,solution.tf));
        end
        subplot(1,nx,n);
        box on
        hold on
        if ~opts.ifexact
            plot(grid,solution.final.primary{2}(:,n),'-d','Color',colors_PDE(n+psize.n0,:),'MarkerSize',marker_size,'LineWidth',line_width,'DisplayName','Numerical solution');
        else
            plot(grid,solution.final.primary{2}(:,n),'d','Color',colors_PDE(n+psize.n0,:),'LineWidth',marker_size,'DisplayName','Numerical solution');
            plot(exact_grid,exsol_grid,'k-','DisplayName','Analytical solution','LineWidth',line_width); 
            legend('FontSize',13,'Interpreter','latex');
        end
        hold off
        xlabel('$s$','FontSize',15,'Interpreter','latex');
        ylabel('$\mathbf{x}$','FontSize',15,'Interpreter','latex');
        if nx==1
            title(['PDE State at Final Time, $\mathbf{x}(t=',num2str(solution.tf),',s)$'],'FontSize',15,'Interpreter','latex');
        else
            title(['$\mathbf{x}_',num2str(n+psize.n0),'(t=',num2str(solution.tf),',s)$'],'FontSize',15,'Interpreter','latex');
        end
        set(gca,'XLim',[min(grid),max(grid)]);
        set(gca,'TickLabelInterpreter','latex');
    end
    figs = [figs,{fig4}];

    if opts.intScheme~=1 | solution.tf==0 | solution.tf==opts.dt
        % Can't plot simulated state evolution in a surf plot
        return
    end

    % % Plot surface plot of state evolution at all points.
    Ngp = size(solution.timedep.primary{2},1);       % Number of grid points
    fig5 = figure('Position',[200 150 fig_width 450]);
    set(gcf, 'Color', 'w');
    box on
    if nx>1
        sgtitle('PDE State Evolution','Interpreter','latex','FontSize',16)
    end
    for ii=1:nx
        x_ii = reshape(solution.timedep.primary{2}(:,ii,t_idcs),Ngp,[]);
        subplot(1,nx,ii);
        surf(tval,grid,x_ii,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
        h = colorbar;
        colormap jet
        box on
        if nx==1
            title('PDE State Evolution $\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
        else
            title(['$\mathbf{x}_',num2str(ii)+psize.n0,'(t,s)$'],'FontSize',15,'Interpreter','latex');
        end
        xlabel('$t$','FontSize',15,'Interpreter','latex');
        ylabel('$s$','FontSize',15,'Interpreter','latex');
        zlabel('$\mathbf{x}$','FontSize',15,'Interpreter','latex');
        set(gca,'XLim',[min(tval),max(tval)]);
        set(gca,'TickLabelInterpreter','latex');
    end
    figs = [figs,{fig5}];

elseif ~strcmp(opts.type,'DDE')
    disp('PDE solution is infinite. Unable to plot. Numerical value is not returned.')
end


end
