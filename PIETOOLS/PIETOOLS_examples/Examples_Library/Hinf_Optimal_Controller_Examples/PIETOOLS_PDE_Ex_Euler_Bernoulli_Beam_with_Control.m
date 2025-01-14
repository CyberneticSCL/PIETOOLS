function PDE_t = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
%
% %---------------------------------------------------------------------% %
% % Euler-Bernoulli Beam Equation with in-domain control
% % PDE        v1_{t}(t,s) = -c*v2_{ss}(t,s) +w(t) +u(t);
% %            v2_{t}(t,s) = v1_{ss}(t,s)
% % With BCs     v1(t,s=0) = v1_{s}(t,s=0) = 0;
% %              v2(t,s=1) = v2_{s}(t,s=1) = 0;
% % Output            z(t) = [u(t); int_{0}^{1}0.5*(1-s)^2*v2(t,s)ds];
% % Parameter c can be set, defaults to 0.1;
% % Example 21 from Shivakumar, 2022 (see bottom of file).
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_control = 1;')

c = 0.1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% pde_var input format
clear stateNameGenerator
v1 = pde_var(s,[0,1]);      v2 = pde_var(s,[0,1]);
w = pde_var('in');          z = pde_var('out',2);
u = pde_var('control');
PDE_t = [diff(v1,'t')==-c*diff(v2,s,2)+w+u;
         diff(v2,'t')==diff(v1,s,2);
         z==[u; int(0.5*(1-s)^2*v2,s,[0,1])];
         subs(v1,s,0)==0;   subs(diff(v1,s),s,0)==0;
         subs(v2,s,1)==0;   subs(diff(v2,s),s,1)==0];


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