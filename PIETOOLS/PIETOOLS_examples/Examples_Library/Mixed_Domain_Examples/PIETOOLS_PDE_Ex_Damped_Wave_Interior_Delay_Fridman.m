function [PDE_t, st] = PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_DAMPED_WAVE_INTERIOR_DELAY_FRIDMAN(GUI,
% PARAMS) returns the damped wave equation with an interior delay below,
% the delay written as a transport state y(s1,s2), as a 'pde_struct'
% object, and (second output) the recommended settings ST of its
% stability test. It sets stability = 1 and stability_form = 'Q'
% (generalized storage T'Q) in the base workspace; a nonzero GUI prints
% that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'tau=1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         (optional) recommended settings (lpisettings_sop
%               structure) for the stability test of this example.
%
% %---------------------------------------------------------------------% %
% % Damped wave equation with an interior delay (Fridman, Orlov 2009,
% % Example 2, with the feedback u = -z(t-tau,s1)):
% % PDE         z_{tt}(t,s1) = 0.1*z_{s1s1} - 2*z_{t} - z(t-tau,s1);
% %                                                       s1 in [0,pi]
% % With BCs    z(t,s1=0) = 0;   z(t,s1=pi) = 0;
% %
% % Let v = z_{t}, and y(t,s1,s2) = z(t-s2,s1), s2 in [0,tau]:
% % PDE         z_{t} = v;
% %             v_{t} = 0.1*z_{s1s1} - 2*v - y(t,s1,s2=tau);
% %             y_{t} = -y_{s2};
% % With BCs    z = y = 0 at s1 = 0, pi;   y(t,s1,s2=0) = z(t,s1).
% % (v carries no derivative in s1.)
% %
% % Parameter tau can be set; default tau = 0.9*tau*.
% %
% % Threshold [derived from the modes]: mode k obeys
% % y'' + 2y' + 0.1*k^2*y + y(t-tau) = 0. A crossing of the imaginary axis
% % requires 0.1*k^2 < 1; the per-mode delay margins are 2.881401 (k = 1),
% % 3.497930 (k = 2), 8.758610 (k = 3), and modes k >= 4 are stable for
% % every delay:
% %   stable iff tau < tau* = 2.881401.
% % Rightmost root of mode 1: -0.1949 at 0.5*tau*, -0.0200 at 0.9*tau*,
% % +0.0155 at 1.1*tau*. The paper reports stability with decay rate 0.8
% % for tau <= 0.31 [read], a different question.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman
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
% Initial coding MMP - 11/17/2026
% MMP, 11/22/2026: st with Q of degree 1 in every direction, slack weights
%   in the directions a space does not own, and the slack integral bases
%   capped at degree Dmax-1: the smallest setting found that certifies
%   stability under an acceptance that also bounds the residual of the
%   solve.

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test

% Default delay (0.9*tau*), overridden by the assignment strings in params
tau = 0.9*2.881401;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% pde_var input
clear stateNameGenerator
z = pde_var(s1,[0,pi]);
v = pde_var(s1,[0,pi],0);                       % v = z_{t}: no derivative in s1
y = pde_var([s1;s2],[0,pi;0,tau],[2;1]);        % 2nd order in s1, 1st in s2
PDE_t = [diff(z,'t')==v;
         diff(v,'t')==0.1*diff(z,s1,2)-2*v-subs(y,s2,tau);
         diff(y,'t')==-diff(y,s2);
         subs([z;y],s1,0)==0;       subs([z;y],s1,pi)==0;
         subs(y,s2,0)==z];


% Example-specific settings (second output only): the light tier, with      % MMP, 11/22/2026
% the changes below                                                         % MMP, 11/22/2026
if nargout>1
    st = lpisettings_sop('light');
    st.lpi.ineq.dD = 0;             % slack lift degree at its minimum Dmin
%   st.sop.deg_Q{2} = [2 1 1];      % Q of degree 2 (multiplier), 1 (integrals) % MMP, 11/22/2026 (was)
    st.sop.deg_Q{2} = [1 1 1];      % Q of degree 1 in every direction      % MMP, 11/22/2026
    st.lpi.ineq.wcross = 1;         % slack weights in unowned directions   % MMP, 11/22/2026
    st.lpi.ineq.cap_inc = 'auto';   % slack integral bases, degree Dmax-1   % MMP, 11/22/2026
end

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{fridman2009exponential,
%   title={Exponential stability of linear distributed parameter systems with time-varying delays},
%   author={Fridman, Emilia and Orlov, Yury},
%   journal={Automatica},
%   volume={45},
%   number={1},
%   pages={194--201},
%   year={2009},
%   doi={10.1016/j.automatica.2008.06.006}
% }
