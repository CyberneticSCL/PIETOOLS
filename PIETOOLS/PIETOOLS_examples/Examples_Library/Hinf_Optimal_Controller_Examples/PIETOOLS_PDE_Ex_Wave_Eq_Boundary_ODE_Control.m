function [PDE_t] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_ODE_Control(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system declared with 
%               pde_var input format.
%
% %---------------------------------------------------------------------% %
% % Wave Equation with boundary control through ODE
% % ODE          xo_{t}(t) = u(t);
% % PDE        x1_{t}(t,s) = x2(t,s);
%              x2_{t}(t,s) = c*x1_{ss}(t,s) + s*(2-s)*w;
% % With BCs     x1(t,s=0) = 0;
% %              x2(t,s=0) = 0;
% %          x1_{s}(t,s=1) = xo(t);
% % and regulated output  z(t) = [xo(t); int(x1(t,s),s,0,1)].
% % Parameter c can be set, defaults to 1.
% % Adapted from Example 23 from Shivakumar, 2022 (see bottom of file).
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s t

%%% Executive Function:
evalin('base','Hinf_control = 1;')

c = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% % % Construct the PDE.
% Standard implementation using phi = [x;x_{t}]
xo = pde_var('state',1,[],[]);
x1 = pde_var(s,[0,1]);          x2 = pde_var(s,[0,1]);
z = pde_var('output',2);        w = pde_var('input',1);
u = pde_var('control',1);
PDE_t = [diff(xo,t)==u;
         diff(x1,t)==x2
         diff(x2,t)==c*diff(x1,s,2)+s*(2-s)*w; 
         subs(x1,s,0)==0; 
         subs(x2,s,0)==0;
         subs(diff(x1,s),s,1)==xo;
         z==[int(x1,s,[0,1]); xo]];

%Alternative implementation using phi = [x;x_{t}-x_{s}]
%  xo = pde_var('state',1,[],[]);
% eta = pde_var('state',1,s,[0,1]);  
% v = pde_var('state',1,s,[0,1]);  
% z = pde_var('output',2);        w = pde_var('input',1);
% u = pde_var('control',1);
% PDE_t=[diff(xo,t)==u;
%        diff(eta,t)==diff(eta,s)+v;
%        diff(v,t)==-diff(v,s)+s*(2-s)*w;
%        z==[xo;int(eta,s,[0,1])];
%        subs(diff(eta,s),s,1)==xo;
%        subs(eta,s,0)==0;
%        subs(v,s,0)==-subs(diff(eta,s),s,0)];

% % Alternative implementation, using phi = [x_{s}; x_{t}]
% phi = pde_var('state',2,s,[0,1]);   x = pde_var('state',1,[],[]);
% w = pde_var('input',1);             r = pde_var('output',1);
% u = pde_var('control');   
% eq_dyn = [diff(x,t,1)==u
%           diff(phi,t,1)==[0 1; c 0]*diff(phi,s,1)+[0;s*(s-1)]*w];
% eq_out= r ==int([1 0]*phi,s,[0,1]);
% bc1 = [0 1]*subs(phi,s,0)==0;   
% bc2 = [1 0]*subs(phi,s,1)==x;
% PDE_t = [eq_dyn;eq_out;bc1;bc2];


if GUI
    disp('No GUI representation available for this system.')
end

end
% @article{shivakumar2022h_,
%   title={$ H\_ $\{$$\backslash$infty$\}$ $-optimal control of coupled ODE-PDE systems using PIE framework and LPIs},
%   author={Shivakumar, Sachin and Das, Amritam and Weiland, Siep and Peet, Matthew},
%   journal={arXiv preprint arXiv:2208.13104},
%   year={2022}
% }