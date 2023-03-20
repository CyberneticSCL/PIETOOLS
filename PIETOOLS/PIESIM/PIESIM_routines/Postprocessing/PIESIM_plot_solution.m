%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_plot_solution.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIESIM_plot_solution(solution, psize, uinput, grid, opts);
% This routine outputs and plots solution of ODE and PDE states


% Inputs: 
% 1) solution 
% solution is a strucutre with the following fields
% --- solution.tf - scalar - actual final time of the solution
% --- solution.final.pde - array of size (N+1) x ns, ns=n0+n1+n2 - pde (distributed state) solution at a final time
% --- solution.final.ode - array of size nx - ode solution at a final time 
% --- solution.final.observed - array of size ny  - final value of observed outputs
% --- solution.final.regulated - array of size nz  - final value of regulated outputs

% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.pde - array of size (N+1) x ns x Nsteps - time-dependent
% --- solution.timedep.ode - array of size nx x Nsteps - time-dependent solution of nx ODE states
%     solution of ns PDE (distributed) states of the primary PDE system
% --- solution.timedep.observed - array of size ny x Nsteps -
%     time-dependent value of observed outputs
% --- solution.timedep.regulated - array of size nz x Nsteps -
%     time-dependent value of regulated outputs

% 2) psize - size of the PIE problem: nw, nu, nf, nx 
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

syms sx;

   % Output solution of ODE states


    if(psize.nx>0)
if (~isempty(solution.final.ode))
for i=1:psize.nx
formatSpec = 'Solution of an ODE state %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.ode(i));
end
else
disp('ODE solution is infinite. Unable to plot. Numerical value is not returned.')
end
    end

% Output values for regulated outputs

    if(psize.nz>0)
if (~isempty(solution.final.regulated))
for i=1:psize.nz
formatSpec = 'Value of a regulated output %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.regulated(i));
end
else
disp('Regulated output is infinite. Unable to plot. Numerical value is not returned.')
end
    end

% Output values for observed outputs

if(psize.ny>0)
if (~isempty(solution.final.observed))
for i=1:psize.ny
formatSpec = 'Value of an observed output %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.observed(i));
end
else
disp('Observed output is infinite. Unable to plot. Numerical value is not returned.')
end
end

% Temporal evolution of solution is only available with BDF scheme
if (opts.intScheme==1 & strcmp(opts.plot,'yes'))
odesol=solution.timedep.ode';
dtime=solution.timedep.dtime;


% Plotting ODE solutions on the same plot
if (psize.nx>0)

    for i=1:psize.nx
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
if (psize.ny>0)
    if (~isempty(solution.timedep.observed))
y=solution.timedep.observed';

    for i=1:psize.ny
     labels{i} = ['y_',num2str(i)];
    end 

 figure;
 plot(dtime,y,'-S','linewidth',2); 
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
    else
disp('Observed output is infinite. Unable to plot. Numerical value is not returned.')
    end
    end

% Plotting regulated outputs on the same plot
if (psize.nz>0)
    if (~isempty(solution.timedep.regulated))
z=solution.timedep.regulated';

    for i=1:psize.nz
     labels{i} = ['z_',num2str(i)];
    end 

 figure;
 plot(dtime,z,'-S','linewidth',2); 
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
else
disp('Regulated output is infinite. Unable to plot. Numerical value is not returned.')
end
end % psize

  if (~strcmp(opts.type,'DDE')& strcmp(opts.plot,'yes'))
    if (~isempty(solution.final.pde))
  if (uinput.ifexact==true)
      a=uinput.a;
      b=uinput.b;
  exact_grid=linspace(a,b,101);
  exsol_grid=double.empty(101,0);
  end % endif (uinput.ifexact)
  ns=sum(psize.n);
  for n=1:ns
      if (uinput.ifexact==true)
  exsol_grid_time=subs(uinput.exact(n),sx,exact_grid);
  exsol_grid=double(subs(exsol_grid_time,solution.tf));
      end % endif (uinput.ifexact)
  figure;
  plot(grid.phys,solution.final.pde(:,n),'rd','MarkerSize',12,'linewidth',2); hold on;
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
    else
disp('PDE solution is infinite. Unable to plot. Numerical value is not returned.')
end

  end % strcmp

end % Temporal evolution


