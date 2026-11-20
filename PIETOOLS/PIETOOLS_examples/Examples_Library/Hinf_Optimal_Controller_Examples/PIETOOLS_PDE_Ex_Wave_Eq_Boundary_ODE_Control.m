function [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_ODE_Control(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_WAVE_EQ_BOUNDARY_ODE_CONTROL(GUI,PARAMS)
% returns the wave equation below, controlled at s = 1 through the ODE
% state xo, with disturbance w and regulated output z, as a 'pde_struct'
% object, and the recommended settings ST of its H-infinity optimal
% controller synthesis. It sets Hinf_control = 1 in the base workspace; a
% nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system declared with 
%               pde_var input format.
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Wave Equation with boundary control through ODE
% % ODE          xo_{t}(t) = u(t);
% % PDE        x1_{t}(t,s) = x2(t,s);
%              x2_{t}(t,s) = c*x1_{ss}(t,s) + s*(2-s)*w;
% % With BCs     x1(t,s=0) = 0;
% %              x2(t,s=0) = 0;
% %          x1_{s}(t,s=1) = xo(t);
% % and regulated output  z(t) = [int(x1(t,s),s,0,1); xo(t)].
% % Parameter c can be set, defaults to 1.
% % Adapted from Example 23 from Shivakumar, 2022 (see bottom of file).
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Wave_Eq_Boundary_ODE_Control
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
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
% function [PDE_t] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_ODE_Control(GUI,params)  % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s t

%%% Executive Function:
evalin('base','Hinf_control = 1;')

% Default parameter, overridden by the assignment strings in params
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


% No GUI representation of this system
if GUI
    disp('No GUI representation available for this system.')
end

% Example-specific settings                                                 % MMP, 11/19/2026
% (veryheavy tier; 1-D storage with upper-integral s-length 2 and
% theta-lengths 3, 4, 4 (multiplier, upper, lower integral); free
% operators of degree 3 in the multiplier and 6 in the integral kernels)
st = lpisettings_sop('veryheavy');                                          % MMP, 11/19/2026
st.lpi.pos.len_com{1} = 2;                                                  % MMP, 11/19/2026
st.lpi.pos.len_sos{1} = [3;4;4];                                            % MMP, 11/19/2026
st.lpi.var.deg{1} = [3 6 6];                                                % MMP, 11/19/2026

end
% @article{shivakumar2022h_,
%   title={$ H\_ $\{$$\backslash$infty$\}$ $-optimal control of coupled ODE-PDE systems using PIE framework and LPIs},
%   author={Shivakumar, Sachin and Das, Amritam and Weiland, Siep and Peet, Matthew},
%   journal={arXiv preprint arXiv:2208.13104},
%   year={2022}
% }