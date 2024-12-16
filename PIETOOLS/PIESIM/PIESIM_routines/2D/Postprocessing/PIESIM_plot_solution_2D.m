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
% it is array of size (N+1) x (nx+ny), nx - number of states dependending on only x, 
%                                      ny - number of states dependending on only y 
% --- solution.final.pde{2} - array containing the solution for states that are the functions of two variables - 
% it is array of size (N+1) x (N+1) x ns, ns - number of states dependending on both x and y
% --- solution.final.ode - array of size no - ode solution at a final time 
% --- solution.final.observed - array of size noo  - final value of observed outputs
% --- solution.final.regulated - array of size nro  - final value of regulated outputs

% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.pde{1,2} - array containing all time-dependent pde solutions 
% --- solution.timedep.pde{1} - array containing the solution for states that are only the functions of one variable - 
% it is array of size (N+1) x (nx+ny) x Nsteps, nx - number of states dependending on only x, 
%                                      ny - number of states dependending on only y 
% --- solution.timedep.pde{2} - array containing the solution for states that are the functions of two variables - 
% it is array of size (N+1) x (N+1) x ns x Nsteps, ns - number of states dependending on both x and y
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
% DJ, 07/17/2024    -   Bugfix for extraction of exact solutions.

syms sx sy;

   % Output solution of ODE states

for i=1:psize.no
formatSpec = 'Solution of an ODE state %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.ode(i));
end

% Output values for regulated outputs
for i=1:psize.nro
formatSpec = 'Value of a regulated output %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.regulated(i));
end

% Output values for observed outputs
for i=1:psize.noo
formatSpec = 'Value of an observed output %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.observed(i));
end

% Temporal evolution of solution is only available with BDF scheme
if (opts.intScheme==1 & opts.tf~=0 & strcmp(opts.plot,'yes'))
odesol=solution.timedep.ode';
dtime=solution.timedep.dtime;


% Plotting ODE solutions on the same plot
if (psize.no>0)

    for i=1:psize.no
     labels{i} = ['State_',num2str(i)];
    end 

 figure;
 if(~isempty(odesol))
 plot(dtime,odesol,'-S','linewidth',2); 
 title('Time evolution of ODE states');
 xlabel('Time');
 ylabel('Value of an ODE state');
 ax = gca;
 ax.FontSize = 16;
 H=gca;
 H.LineWidth=3;
 [leg, hobj]=legend(labels);
 set(leg,'Box','off','Location','north','FontSize',16);
 hl = findobj(hobj,'type','line');
 set(hl,'LineWidth',2);  
 clear labels;
 else
     disp('ODE solution is infinite. Unable to plot. Numerical value is not returned.')
 end
end

% Plotting observed outputs on the same plot
if (psize.noo>0)
y=solution.timedep.observed';

    for i=1:psize.noo
     labels{i} = ['y_',num2str(i)];
    end 

 figure;
 plot(dtime,y,'-S'); 
 title('Time evolution of observed outputs');
 xlabel('Time');
 ylabel('Value of an observed output');
 ax = gca;
 ax.FontSize = 16;
 H=gca;
 H.LineWidth=3;
 [leg, hobj]=legend(labels);
 set(leg,'Box','off','Location','north','FontSize',16);
 hl = findobj(hobj,'type','line');
 set(hl,'LineWidth',2);  
 clear labels;
end

% Plotting regulated outputs on the same plot
if (psize.nro>0)
z=solution.timedep.regulated';

    for i=1:psize.nro
     labels{i} = ['z_',num2str(i)];
    end 

 figure;
 plot(dtime,z,'-S'); 
 title('Time evolution of regulated outputs');
 xlabel('Time');
 ylabel('Value of a regulated output');
 ax = gca;
 ax.FontSize = 16;
 H=gca;
 H.LineWidth=3;
 [leg, hobj]=legend(labels);
 set(leg,'Box','off','Location','north','FontSize',16);
 hl = findobj(hobj,'type','line');
 set(hl,'LineWidth',2);  
 clear labels;
end
end

% 1D states, differentiable in x
ns_tot = psize.no;
if (sum(psize.nx)>0)
    if (~strcmp(opts.type,'DDE') && strcmp(opts.plot,'yes'))
        if (~isempty(solution.final.pde))
            if (uinput.ifexact==true)
                a=uinput.dom(1,1);
                b=uinput.dom(1,2);
                exact_grid=linspace(a,b,101);
                exsol_grid=double.empty(101,0);
            end % endif (uinput.ifexact)
            ns=sum(psize.nx);
            for n=1:ns
                % Plot the exact solution if desired.
                if (uinput.ifexact==true)
                    exsol_grid_time=subs(uinput.exact(n+ns_tot),sx,exact_grid);
                    exsol_grid=double(subs(exsol_grid_time,solution.tf));
                end % endif (uinput.ifexact)
                figure;
                plot(grid.phys(:,1),solution.final.pde{1}(:,n),'rd','MarkerSize',12,'linewidth',2); hold on;
                if (uinput.ifexact==true)
                plot(exact_grid,exsol_grid,'k','linewidth',3); 
                end % endif (uinput.ifexact)
                xlabel('Spatial variable');
                ylabel('Solution at a final time');
                title('Plot of a primary state solution',num2str(n));
                
                if (uinput.ifexact==true)
                legend('Numerical solution','Analytical solution');
                else
                legend('Numerical solution');
                end % endif (uinput.ifexact)
                ax = gca;
                ax.FontSize = 24;
                H=gca;
                H.LineWidth=3;
            end % for n=1:ns
            % Keep track of total number of state variables
            ns_tot = ns_tot + ns;
        else
            disp('PDE solution is infinite. Unable to plot. Numerical value is not returned.')
        end
    end
end


% 1D states, differentiable in y
if (sum(psize.ny)>0)
    if (~strcmp(opts.type,'DDE') && strcmp(opts.plot,'yes'))
        if (~isempty(solution.final.pde))
            if (uinput.ifexact==true) 
                c=uinput.dom(2,1);
                d=uinput.dom(2,2);
                exact_grid=linspace(c,d,101);
                exsol_grid=double.empty(101,0);
            end % endif (uinput.ifexact)
            ns=sum(psize.ny);
            for n=1:ns
                % Plot the exact solution if desired.
                if (uinput.ifexact==true)
                    exsol_grid_time=subs(uinput.exact(ns_tot+n),sy,exact_grid);
                    exsol_grid=double(subs(exsol_grid_time,solution.tf));
                end % endif (uinput.ifexact)
            figure;
            ng=sum(psize.nx)+n;
            plot(grid.phys(:,2),solution.final.pde{1}(:,ng),'rd','MarkerSize',12,'linewidth',2); hold on;
            if (uinput.ifexact==true)
            plot(exact_grid,exsol_grid,'k','linewidth',3); 
            end % endif (uinput.ifexact)
            xlabel('Spatial variable');
            ylabel('Solution at a final time');
            title('Plot of a primary state solution',num2str(ng));
            
            if (uinput.ifexact==true)
            legend('Numerical solution','Analytical solution');
            else
            legend('Numerical solution');
            end % endif (uinput.ifexact)
            ax = gca;
            ax.FontSize = 24;
            H=gca;
            H.LineWidth=3;
            end % for n=1:ns
            % Keep track of total number of state variables
            ns_tot = ns_tot + ns;
        else
            disp('PDE solution is infinite. Unable to plot. Numerical value is not returned.')
        end
    end
end


% 2D states
if (~strcmp(opts.type,'DDE') && strcmp(opts.plot,'yes'))
    if (uinput.ifexact==true)
    a=uinput.dom(1,1);
    b=uinput.dom(1,2);
    c=uinput.dom(2,1);
    d=uinput.dom(2,2);
    exact_grid_x=linspace(a,b,101);
    exact_grid_y=linspace(c,d,101);
    exsol_grid=double.empty(101,0);
    end
    ns=sum(psize.n);
    for n=1:ns
    if (uinput.ifexact==true)
    exsol_grid_time=subs(subs(uinput.exact(n+ns_tot),sx,exact_grid_x'),sy,exact_grid_y);
    exsol_grid=double(subs(exsol_grid_time,solution.tf));
    exsol_numgrid_time=subs(subs(uinput.exact(n+ns_tot),sx,grid.phys(:,1)),sy,grid.phys(:,2)');
    exsol_numgrid=double(subs(exsol_numgrid_time,solution.tf));
    end
    
    
    % Plot isosurface for numerical solution

    % PIESIM stores nD solutions by columns, i.e., e.g., for u(x,y), the
    % entries with fixed y and varying x are stored in columnes, and the
    % entries with fixed x and varying y are stored in rows

    % MATLAb surf routine does the opposite, so we plot the transpose of
    % the solution vector to get a correct orienation in space 

    % Solution is stored
    figure;
    surf(grid.phys(:,1),grid.phys(:,2),solution.final.pde{2}(:,:,n)');
    xlabel('x'), ylabel('y'), zlabel('Numerical Solution');
    ax = gca;
    ax.FontSize = 24;
    H=gca;
    H.LineWidth=3;
    ng=sum(psize.nx)+sum(psize.ny)+n;
    title('Plot of primary state solution',num2str(ng));

    if (uinput.ifexact==true)
    % Plot isosurface for analytical solution
    figure;
    surf(exact_grid_x,exact_grid_y,exsol_grid');
    xlabel('x'), ylabel('y'), zlabel('Analytical Solution');
    ax = gca;
    ax.FontSize = 24;
    H=gca;
    H.LineWidth=3;
    title('Exact primary state solution',num2str(ng));
    
    figure;
    % Plot isosurface for the diffefence between analytical and numerical solution
    surf(grid.phys(:,1),grid.phys(:,2),exsol_numgrid'-solution.final.pde{2}(:,:,n)');
    xlabel('x'), ylabel('y'), zlabel('Error');
    title('Error in primary state solution',num2str(ng));
    end
    
    ax = gca;
    ax.FontSize = 24;
    H=gca;
    H.LineWidth=3;
    end
end