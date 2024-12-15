function figs = PIESIM_plotsolution(solution,grids,varargin)
% FIGS = PIESIM_PLOTSOLUTION(SOLUTION,GRIDS,VARARGIN) takes a solution from
% PIESIM along with information on the grid on which simulation was
% performed, and plots the simulated evolution of the ODE state components, 
% PDE state components, and regulated outputs against time, returning a
% cell with figure objects.
%
% INPUTS
% - solution:   struct corresponding to the first output of a call to
%               PIESIM, specifying simulation results for some ODE-PDE or
%               PIE system;
% - grids:      struct corresponding to the second output of a call to
%               PIESIM, specifying information on the grid on which
%               simulation was performed. In particular, the spatial
%               positions of the grid points must be specified under
%               "grids.phys";
% - varargin:   Optional additional arguments specifying options for the
%               plot. Each additional argument must be one of the following:
%       'tlog': To plot the evolution against time on log scale;
%   'use_grid': To plot the PDE state evolution at given grid points,
%               rather than plotting a surface plot of the evolution at all
%               positions. If this argument is passed, the subsequent
%               argument must be an array of integers, speicfying the
%               indices of the grid points at which to plot the solution.
%       'name': To add a name to the title of the plots. The subsequent
%               argument must be an element of type char, which will then
%               be added at the start of the titles of each of the returned
%               plots.
%
% OUTPUTS
% - figs:       1x3 cell array, with each element an object of type
%               'Figure'. The first figure gives a plot of the ODE state
%               dynamics, the second of the PDE state dynamics, and the
%               third of the regulated outputs. If the system involves no
%               ODE states, the first element will be excluded. If the
%               system involves no regulated outputs, the last element will
%               be excluded.
%
% NOTES
% The function is currently only supported for 1D simulations.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_plotsolution
%
% Copyright (C)2024  PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 12/15/2024: Initial coding;


% % Check the optional input arguments
use_tlog = false;
use_grid = false;
plot_name = '';
skip_ii = false;
for ii=1:numel(varargin)
    if skip_ii
        skip_ii = false;
        continue
    end
    if strcmpi(varargin{ii},'tlog')
        % Plot time on log scale;
        use_tlog = true;
    elseif strcmpi(varargin{ii},'gridpoints')
        % Plot PDE states only at certain grid points.
        use_grid = true;
        if numel(varargin)>ii && isnumeric(varargin{ii+1})
            grid_idcs = varargin{ii+1};
            skip_ii = true;
        else
            grid_idcs = [];
        end
    elseif strcmpi(varargin{ii},'title')
        % Add a name to the plots at the start of the title.
        if numel(varargin)==ii || (~ischar(varargin{ii+1}) && ~isstring(varargin{ii+1}))
            error('A title has not been properly specified for the plots.')
        else
            plot_name = [char(varargin{ii+1}),' '];
        end
        skip_ii = true;
    end
end

% Determine time indices at which to plot solution.
tval = solution.timedep.dtime;
Nplot = min(length(tval),100);  % plot at at most 100 time points.
if use_tlog
    t_idcs = unique(floor(logspace(0,log(length(tval))/log(10),Nplot)));
    t_idcs = t_idcs(2:end);     % Can't plot at t=0;
else
    t_idcs = floor(linspace(1,length(tval),Nplot));
end
tval = tval(t_idcs);

% Initialize empty output cell.
figs = {};


% % First, plot any ODE state solutions.
n_ODE = size(solution.timedep.ode,1);
if n_ODE>0
    fig1 = figure('Position',[200 150 800 450]);
    set(gcf, 'Color', 'w');
    box on
    title([plot_name,'ODE State Evolution'],'FontSize',16,'Interpreter','latex');
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$X$','FontSize',15,'Interpreter','latex');
    hold on
    for ii=1:n_ODE
        X_ii = solution.timedep.ode(ii,t_idcs);
        plot(tval,X_ii,'DisplayName',['$X_',num2str(ii),'(t)$'],'LineWidth',2);
    end
    hold off
    legend('FontSize',15,'Interpreter','latex');
    if use_tlog
        set(gca,'XScale','log');
    end
    set(gca,'XLim',[min(tval),max(tval)]);
    figs = [figs,{fig1}];
end


% % Next, plot the PDE state solutions.
colors_PDE = [0 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.9290 0.6940 0.1250;
            0.4940 0.1840 0.5560;
            0.4660 0.6740 0.1880;
            0.3010 0.7450 0.9330;
            0.6350 0.0780 0.1840];
n_PDE = size(solution.timedep.pde,2);
if n_PDE>0
    N = size(solution.timedep.pde,1);       % Number of grid points
    fig_width = min(n_PDE*600,1800);
    fig2 = figure('Position',[200 150 fig_width 450]);
    set(gcf, 'Color', 'w');
    box on
    if n_PDE>1
        sgtitle([plot_name,'PDE State Evolution'],'Interpreter','latex','FontSize',16)
    end
    if use_grid
        % % Plot state evolution only at certain grid points.
        if isempty(grid_idcs)
            grid_idcs = (1:N);
        end
        grid_pts = grids.phys(grid_idcs);
        for ii=1:n_PDE
            x_ii = reshape(solution.timedep.pde(:,ii,t_idcs),N,[]);
            ax_ii = subplot(1,n_PDE,ii);
            box on
            hold on
            for jj=1:length(grid_idcs)
                % Plot state component ii at grid point jj;
                x_jj = x_ii(jj,:);
                if jj<=size(colors_PDE,1)-n_ODE
                    plot(tval,x_jj,'Color',colors_PDE(jj+n_ODE,:),'DisplayName',['$s=',num2str(grid_pts(jj)),'$'],'LineWidth',2);
                else
                    plot(tval,x_jj,'DisplayName',['$s=',num2str(grid_pts(jj)),'$'],'LineWidth',2);
                end
            end
            hold off
            title(['$x_',num2str(ii),'(t,s)$'],'FontSize',15,'Interpreter','latex');
            xlabel('$t$','FontSize',15,'Interpreter','latex');
            ylabel(['$x_',num2str(ii),'$'],'FontSize',15,'Interpreter','latex');
            legend('FontSize',13,'Interpreter','latex');
            if use_tlog
                ax_ii.XScale = 'log';
            end
            set(gca,'XLim',[min(tval),max(tval)]);
        end
    else
        % % Plot surface plot of state at all points.
        for ii=1:n_PDE
            x_ii = reshape(solution.timedep.pde(:,ii,t_idcs),N,[]);
            ax_ii = subplot(1,n_PDE,ii);
            surf(tval,grids.phys,x_ii,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
            h = colorbar;
            colormap jet
            box on
            if n_PDE==1
                title([plot_name,'PDE State Evolution $x_',num2str(ii),'(t,s)$'],'FontSize',15,'Interpreter','latex');
            else
                title(['$x_',num2str(ii),'(t,s)$'],'FontSize',15,'Interpreter','latex');
            end
            xlabel('$t$','FontSize',15,'Interpreter','latex');
            ylabel('$s$','FontSize',15,'Interpreter','latex');
            zlabel(['$x_',num2str(ii),'$'],'FontSize',15,'Interpreter','latex');
            if use_tlog
                ax_ii.XScale = 'log';
            end
            set(gca,'XLim',[min(tval),max(tval)]);
        end
    end
    figs = [figs,{fig2}];
end


% % Finally, plot any regulated outputs
n_z = size(solution.timedep.regulated,1);
colors_z = {'b','g','m','r','k','c','r','y'};
if n_z>0
    fig3 = figure('Position',[200 150 800 450]);
    box on
    title([plot_name,'Regulated Output Evolution'],'FontSize',16,'Interpreter','latex');
    xlabel('$t$','FontSize',15,'Interpreter','latex');
    ylabel('$z$','FontSize',15,'Interpreter','latex');
    hold on
    for ii=1:n_z
        z_ii = solution.timedep.regulated(ii,t_idcs);
        if ii<=length(colors_z)
            plot(tval,z_ii,[colors_z{ii},'-'],'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',2);   
        else
            plot(tval,z_ii,'DisplayName',['$z_',num2str(ii),'(t)$'],'LineWidth',2);
        end
    end
    hold off
    legend('FontSize',15,'Interpreter','latex');
    if use_tlog
        set(gca,'XScale','log');
    end
    set(gca,'XLim',[min(tval),max(tval)]);
    set(gcf, 'Color', 'w');
    figs = [figs,{fig3}];
end

end