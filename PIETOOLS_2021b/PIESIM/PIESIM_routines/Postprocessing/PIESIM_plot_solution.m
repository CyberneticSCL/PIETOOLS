%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_plot_solution.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIESIM_plot_solution(solution, psize, uinput, grid, opts);
% This routine outputs and plots solution of ODE and PDE states


% Inputs: 
% 1) solution
% This field contains
% AVAILABLE FOR ALL OPTS.INTSCHEME OPTIONS
% solution.tf - actual final time of the solution
% solution.final.pde - pde (distributed state) solution at a final time : matrix of
% dimension (N+1) x ns, ns=n0+n1+n2
% solution.final.ode - ode solution at a final time : array of dimensions
% nx x 1

% AVAILABLE ONLY FOR OPTS.INTSCHEME=1 (BDF) OPTION
% solution.timedep.ode - array of size nx x Nsteps - time-dependent solution of nx ODE states
% solution.timedep.pde - array of size (N+1) x ns x Nsteps - time-dependent
% solution of ns PDE (distribited) states of the primary PDE system
% solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution

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

syms sx;

   % Output solution of ODE states

for i=1:psize.nx
formatSpec = 'Solution of an ODE state %s at a final time %f is %8.4f\n';
fprintf(formatSpec,num2str(i),solution.tf, solution.final.ode(i));
end

% Temporal evolution of solution is only available with BDF scheme
if (opts.intScheme==1 & strcmp(opts.plot,'yes'))
odesol=solution.timedep.ode';
dtime=solution.timedep.dtime;


% Plotting ODE solutions on the same plot
if (psize.nx>0)

    for i=1:psize.nx
     labels{i} = ['S_',num2str(i)];
    end 

 figure;
 plot(dtime,odesol,'-S'); 
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

end
end

  if (~strcmp(opts.type,'DDE'))
  if (uinput.ifexact==true)
      a=uinput.a;
      b=uinput.b;
  exact_grid=linspace(a,b,101);
  exsol_grid=double.empty(101,0);
  end
  ns=sum(psize.n);
  for n=1:ns
      if (uinput.ifexact==true)
  exsol_grid_time=subs(uinput.exact(n),sx,exact_grid);
  exsol_grid=double(subs(exsol_grid_time,solution.tf));
      end
  figure;
  plot(grid.phys,solution.final.pde(:,n),'rd','MarkerSize',12,'linewidth',2); hold on;
  if (uinput.ifexact==true)
  plot(exact_grid,exsol_grid,'k','linewidth',3); 
  end
  xlabel('Spatial variable');
  ylabel('Solution at a final time');
  title('Plot of a primary state solution',num2str(n));
  ax = gca;
  ax.FontSize = 16;
  H=gca;
  H.LineWidth=3;
  end
  end


