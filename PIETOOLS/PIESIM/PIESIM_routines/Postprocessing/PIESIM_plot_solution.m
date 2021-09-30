%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_plot_solution.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIESIM_plot_solution(solution, psize, uinput, grid, opts);
syms sx;
% This routine outputs and plots solution of ODE and PDE states

% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_19_2021

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
 figure;
 all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'}; 
 all_colors={'red', 'green', 'blue', 'cyan', 'magenta', 'yellow','black'};
for i=1:psize.nx
  plot(dtime(:),odesol(:,i),'Marker',all_marks{mod(i,13)},'Color',all_colors{mod(i,7)}, 'MarkerSize',6,'linewidth',1); hold on;  
  formatSpec = 'ODE state %s has marker %s and color %s\n';
  fprintf(formatSpec, num2str(i), all_marks{mod(i,13)}, all_colors{mod(i,7)});
  xlabel('Time');
  ylabel('Value of an ODE state');
  title('Plot of ODE state solution(s)');
  ax = gca;
  ax.FontSize = 16;
  H=gca;
  H.LineWidth=3;
  hold on;
end
  

% Plotting ODE solutions on different plots
% for i=1:psize.nx
%  figure;
%   plot(dtime(:),odesol(:,i),'rd','MarkerSize',6,'linewidth',1); hold on; 
%   xlabel('Time');
%   ylabel('Value of an ODE state');
%   title('Plot of an ODE state solution',num2str(i));
%   ax = gca;
%   ax.FontSize = 16;
%   H=gca;
%   H.LineWidth=3;
%   hold on;
% end
end
end


 
  if (~strcmp(opts.type,'DDE')&strcmp(opts.plot,'yes'))
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
  
  
