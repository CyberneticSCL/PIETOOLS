%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_plot_solution.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figs = PIESIM_plot_solution(solution, psize, uinput, grid, opts);
% This routine outputs and plots solution of ODE and PDE states


% Inputs: 
% 1) solution 
% solution is a strucutre with the following fields
% --- solution.tf - scalar - actual final time of the solution
% --- solution.final.pde - array of size (N+1) x ns, ns=n0+n1+n2 - pde (distributed state) solution at a final time
% --- solution.final.ode - array of size no - ode solution at a final time 
% --- solution.final.observed - array of size noo  - final value of observed outputs
% --- solution.final.regulated - array of size nro  - final value of regulated outputs

% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.pde - array of size (N+1) x ns x Nsteps - time-dependent
% --- solution.timedep.ode - array of size no x Nsteps - time-dependent solution of no ODE states
% --- solution.timedep.observed - array of size noo x Nsteps -
%     time-dependent value of observed outputs
% --- solution.timedep.regulated - array of size nro x Nsteps -
%     time-dependent value of regulated outputs

% 2) psize - size of the PIE problem: nw, nu, nf, no 
% 3) uinput -  user-defined boundary inputs, forcing functions and initial conditions
% 4) grid - physical and computational grid for n0 states
% 5) opts - options for simulation parameters
% 
% Outputs: none


% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_19_2021
% YP 6/16/2022 - added outputs for observed and regulated outputs 
% DJ, 12/26/2024 - Plot PDE state evolution using surf;
% DJ, 01/01/2025 - Change plot specifications;

syms sx;

   
% % Display final value of ODE state, if possible.
if(psize.no>0)
    if (~isempty(solution.final.ode))
        for i=1:psize.no
            formatSpec = 'Solution of ODE state %s at a final time %f is %8.4f\n';
            fprintf(formatSpec,num2str(i),solution.tf, solution.final.ode(i));
        end
        plot_ode = true;
    else
        disp('ODE solution is infinite. Unable to plot. Numerical value is not returned.')
        plot_ode = false;
    end
end

% % Display final value of regulated output, if possible.
if(psize.nro>0)
    if (~isempty(solution.final.regulated))
        for i=1:psize.nro
            formatSpec = 'Value of regulated output %s at a final time %f is %8.4f\n';
            fprintf(formatSpec,num2str(i),solution.tf, solution.final.regulated(i));
        end
        plot_z = true;
    else
        disp('Regulated output is infinite. Unable to plot. Numerical value is not returned.')
        plot_z = false;
    end
end

% % Display final value of observed output, if possible.
if(psize.noo>0)
    if (~isempty(solution.final.observed))
        for i=1:psize.noo
            formatSpec = 'Value of observed output %s at a final time %f is %8.4f\n';
            fprintf(formatSpec,num2str(i),solution.tf, solution.final.observed(i));
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
    disp("Final time is t=0; no plots produced.")
    return
else
    % Determine time indices at which to plot solution.
    dtime = solution.timedep.dtime;
    Nplot = min(length(dtime),100);  % plot at at most 100 time points.
    t_idcs = floor(linspace(1,length(dtime),Nplot));
    tval = dtime(t_idcs);
end
figs = {};


% % Plot the ODE states evolution in a single plot
if (psize.no>0) && plot_ode && opts.intScheme==1
    odesol = solution.timedep.ode';

    fig1 = figure('Position',[200 150 800 450]);
    set(gcf, 'Color', 'w');
    box on
    
    hold on
    for ii=1:psize.no
        X_ii = odesol(t_idcs,ii)';
        plot(tval,X_ii,'DisplayName',['$x_',num2str(ii),'(t)$'],'LineWidth',2);
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
    if psize.no>1
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
if (psize.noo>0) && plot_y && opts.intScheme==1

    fig2 = figure('Position',[200 150 800 450]);
    box on
    title('Observed Output Evolution','FontSize',16,'Interpreter','latex');
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$y$','FontSize',15,'Interpreter','latex');
    hold on
    for ii=1:psize.noo
        y_ii = solution.timedep.observed(ii,t_idcs);
        if ii<=size(colors_y,1)
            plot(tval,y_ii,'Color',colors_y(psize.noo-ii+1,:),'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',2);   
        else
            plot(tval,y_ii,'DisplayName',['$y_',num2str(ii),'(t)$'],'LineWidth',2);
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
end

% % Plot the regulated outputs evolution in a single plot.
colors_z = {'b','g','m','r','k','c','r','y'};
if (psize.nro>0) && plot_z && opts.intScheme==1

    fig3 = figure('Position',[200 150 800 450]);
    box on
    title('Regulated Output Evolution','FontSize',16,'Interpreter','latex');
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$z$','FontSize',15,'Interpreter','latex');
    hold on
    for ii=1:psize.nro
        z_ii = solution.timedep.regulated(ii,t_idcs);
        if ii<=length(colors_z)
            plot(tval,z_ii,[colors_z{ii},'-'],'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',2);   
        else
            plot(tval,z_ii,'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',2);
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
end

% % Plot the PDE states
if ~strcmp(opts.type,'DDE') && ~isempty(solution.final.pde)

    colors_PDE = [0 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.9290 0.6940 0.1250;
            0.4940 0.1840 0.5560;
            0.4660 0.6740 0.1880;
            0.3010 0.7450 0.9330;
            0.6350 0.0780 0.1840];

    % Initialize grid points at which to compute exact solution.
    if (uinput.ifexact==true)
        a = uinput.a;
        b = uinput.b;
        exact_grid = linspace(a,b,round(Nplot*3/4));
        exsol_grid = double.empty(round(Nplot*3/4),0);
    end
    ns = sum(psize.n);

    % Plot the PDE states at the final time
    fig_width = min(ns*600,1800);
    fig4 = figure('Position',[200 150 fig_width 450]);
    set(gcf, 'Color', 'w');
    box on
    if ns>1
        sgtitle('PDE State at Final Time','Interpreter','latex','FontSize',16)
    end
    for n=1:ns
        if (uinput.ifexact==true)
            exsol_grid_time = subs(uinput.exact(n),sx,exact_grid);
            exsol_grid = double(subs(exsol_grid_time,solution.tf));
        end
        subplot(1,ns,n);
        box on
        hold on
        if ~uinput.ifexact
            plot(grid.phys,solution.final.pde(:,n),'-d','Color',colors_PDE(n+psize.no,:),'LineWidth',2,'DisplayName','Numerical solution');
        else
            plot(grid.phys,solution.final.pde(:,n),'d','Color',colors_PDE(n+psize.no,:),'LineWidth',2,'DisplayName','Numerical solution');
            plot(exact_grid,exsol_grid,'k-','DisplayName','Analytical solution'); 
            legend('FontSize',13,'Interpreter','latex');
        end
        hold off
        xlabel('$s$','FontSize',15,'Interpreter','latex');
        ylabel('$\mathbf{x}$','FontSize',15,'Interpreter','latex');
        if ns==1
            title(['PDE State at Final Time, $\mathbf{x}(t=',num2str(solution.timedep.dtime(end)),',s)$'],'FontSize',15,'Interpreter','latex');
        else
            title(['$\mathbf{x}_',num2str(n+psize.no),'(t=',num2str(solution.timedep.dtime(end)),',s)$'],'FontSize',15,'Interpreter','latex');
        end
        set(gca,'XLim',[min(grid.phys),max(grid.phys)]);
        set(gca,'TickLabelInterpreter','latex');
    end
    figs = [figs,{fig4}];

    if opts.intScheme~=1
        % Can't plot simulated state evolution
        return
    end

    % % Plot surface plot of state evolution at all points.
    N = size(solution.timedep.pde,1);       % Number of grid points
    fig5 = figure('Position',[200 150 fig_width 450]);
    set(gcf, 'Color', 'w');
    box on
    if ns>1
        sgtitle('PDE State Evolution','Interpreter','latex','FontSize',16)
    end
    for ii=1:ns
        x_ii = reshape(solution.timedep.pde(:,ii,t_idcs),N,[]);
        subplot(1,ns,ii);
        surf(tval,grid.phys,x_ii,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
        h = colorbar;
        colormap jet
        box on
        if ns==1
            title('PDE State Evolution $\mathbf{x}(t,s)$','FontSize',15,'Interpreter','latex');
        else
            title(['$\mathbf{x}_',num2str(ii)+psize.no,'(t,s)$'],'FontSize',15,'Interpreter','latex');
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