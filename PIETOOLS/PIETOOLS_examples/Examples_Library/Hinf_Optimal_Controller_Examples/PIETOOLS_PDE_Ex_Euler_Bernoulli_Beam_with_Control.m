function [PDE_t,st] = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_EULER_BERNOULLI_BEAM_WITH_CONTROL(GUI,PARAMS)
% returns the Euler-Bernoulli beam below, with distributed control u,
% disturbance w and regulated output z, as a 'pde_struct' object, and the
% recommended settings ST of its H-infinity optimal controller synthesis.
% It sets Hinf_control = 1 in the base workspace; a nonzero GUI prints
% that no GUI representation exists.
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
% % Euler-Bernoulli Beam Equation with in-domain control
% % PDE        v1_{t}(t,s) = -c*v2_{ss}(t,s) +s^2*w(t) +u(t,s);
% %            v2_{t}(t,s) = v1_{ss}(t,s)
% % With BCs     v1(t,s=0) = v1_{s}(t,s=0) = 0;
% %              v2(t,s=1) = v2_{s}(t,s=1) = 0;
% % Output            z(t) = int_{0}^{1} (1-s)*v2(t,s) ds;
% % Parameter c can be set, defaults to 0.1;
% % Adapted from Example 21 from Shivakumar, 2022 (see bottom of file).
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control
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
%   help revised (regulated output z corrected to the one the code sets),
%   license added.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function PDE_t = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_control = 1;')

% Default parameter, overridden by the assignment strings in params
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
u = pde_var('control',s,[0,1]);
w = pde_var('in');          z = pde_var('out');
PDE_t = [diff(v1,'t')==-c*diff(v2,s,2)+u+s^2*w;
         diff(v2,'t')==diff(v1,s,2);
         z==int((1-s)*v2,s,[0,1]);
         subs(v1,s,0)==0;   subs(diff(v1,s),s,0)==0;
         subs(v2,s,1)==0;   subs(diff(v2,s),s,1)==0];


% No GUI representation of this system
if GUI
    disp('No GUI representation available for this system.')
end

% Example-specific settings                                                 % MMP, 11/21/2026
% (veryheavy tier; 1-D storage with s-lengths 2 (lower integral), 3 (upper
% integral) and theta-lengths 4, 4, 3 (multiplier, upper, lower integral);
% slack lift degree at its minimum Dmin; free operators of degree 1 in the
% multiplier and 6, 2 in the integral kernels (output, integration
% variable); no slack weight raise of the synthesis LPI)
st = lpisettings('veryheavy');                                              % MMP, 11/21/2026
st.lpi.pos.len_inc{1} = 2;                                                  % MMP, 11/21/2026
st.lpi.pos.len_com{1} = 3;                                                  % MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [4;4;3];                                            % MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/21/2026
st.lpi.var.deg{1} = [1 2 6];                                                % MMP, 11/21/2026
st.sop.dw_P = 0;                                                            % MMP, 11/21/2026

end
% @article{shivakumar2022h_,
%   title={$ H\_ $\{$$\backslash$infty$\}$ $-optimal control of coupled ODE-PDE systems using PIE framework and LPIs},
%   author={Shivakumar, Sachin and Das, Amritam and Weiland, Siep and Peet, Matthew},
%   journal={arXiv preprint arXiv:2208.13104},
%   year={2022}
% }