function [PDE_t,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Boundary_ODE_Control(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_REACTION_DIFFUSION_BOUNDARY_ODE_CONTROL(GUI,
% PARAMS) returns the reaction-diffusion equation below, controlled at
% s = 1 through the ODE state x, with disturbance w and regulated output
% z, as a 'pde_struct' object, and the recommended settings ST of its
% H-infinity optimal controller synthesis. It sets Hinf_control = 1 in the
% base workspace; a nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
% - st:         example-specific settings (lpisettings structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Example reaction-diffusion equation with boundary control through ODE
% % ODE         x_{t}(t) = u(t);
% % PDE       v_{t}(t,s) = lam*v(t,s) + v_{ss}(t,s) + s*(2-s)*w(t)
% % With BCs    v(t,s=0) = 0
% %             v_{s}(t,s=1) = x(t)
% % Output      z(t) = [x(t); int_{0}^{1}v(t,s)ds]
% % Parameter lam can be set, defaults to 5.
% % Optimal state feedback can achieve L2 gain bound of 4.99
% % Adapted from Example 22 from Shivakumar, 2022 (see bottom of file).
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Reaction_Diffusion_Boundary_ODE_Control
%
% Copyright (C) 2026 PIETOOLS Team
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DJ - 01/13/2025
% MMP, 11/21/2026: Changes since PIETOOLS 2025: second output st, the
%   recommended settings of the H-infinity optimal controller synthesis;
%   help revised, license added.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Boundary_ODE_Control(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_control = 1;')

% Default parameter, overridden by the assignment strings in params
lam = 5;
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
pde_var state x v  input w  output z1 z2  control u
v.vars = s;     v.dom = [0,1];
PDE_t = [diff(x,'t')==u;
         diff(v,'t')==lam*v+diff(v,s,2)+s*(2-s)*w;
         z1==x;
         z2==int(v,s,[0,1]);
         subs(v,s,0)==0;
         subs(diff(v,s),s,1)==x];


% No GUI representation of this system
if GUI
    disp('No GUI representation available for this system.')
end

% Example-specific settings                                                 % MMP, 11/21/2026
% (light tier; 1-D storage with upper-integral s-length 1 and theta-lengths
% 3, 2, 1 (multiplier, upper, lower integral); slack integral bases capped
% at total degree theta+s <= Dmax-1; free operators of degree 0 in the
% multiplier and 2 in the integral kernels)
st = lpisettings('light');                                                  % MMP, 11/21/2026
st.lpi.pos.len_com{1} = 1;                                                  % MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [3;2;1];                                            % MMP, 11/21/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/21/2026
st.lpi.var.deg{1} = [0 2 2];                                                % MMP, 11/21/2026

end
% @article{shivakumar2022h_,
%   title={$ H\_ $\{$$\backslash$infty$\}$ $-optimal control of coupled ODE-PDE systems using PIE framework and LPIs},
%   author={Shivakumar, Sachin and Das, Amritam and Weiland, Siep and Peet, Matthew},
%   journal={arXiv preprint arXiv:2208.13104},
%   year={2022}
% }