function [PDE_t, st] = PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_DAMPED_WAVE_INTERIOR_DELAY_RIEMANN_FRIDMAN(
% GUI,PARAMS) returns the damped wave equation with an interior delay
% below in first-order (Riemann) form on the half interval, the delay
% written as a transport state y(s1,s2), as a 'pde_struct' object, and
% (second output) the recommended settings ST of its stability test. It
% sets stability = 1 and stability_form = 'Q' (generalized storage T'Q)
% in the base workspace; a nonzero GUI prints that no GUI representation
% exists.
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
% - st:         (optional) recommended settings (lpisettings
%               structure) for the stability test of this example.
%
% %---------------------------------------------------------------------% %
% % Damped wave equation with an interior delay (Fridman, Orlov 2009,
% % Example 2), first-order (Riemann) form on the half interval:
% % PDE         z_{tt}(t,s) = 0.1*z_{ss} - 2*z_{t} - z(t-tau,s);
% %                                                       s in [0,pi]
% % With BCs    z(t,s=0) = 0;   z(t,s=pi) = 0;
% %
% % Restricted to the modes symmetric about s = pi/2 (the odd sine modes,
% % which contain the critical mode k = 1): s1 in [0,pi/2], z(t,0) = 0,
% % z_{s}(t,pi/2) = 0. Let e = z_{s}, v = z_{t}, and
% % y(t,s1,s2) = e(t-s2,s1), s2 in [0,tau], so that
% % z(t-tau,s1) = int_0^s1 y(t,r,tau) dr:
% % PDE         e_{t} = v_{s1};
% %             v_{t} = 0.1*e_{s1} - 2*v - int_0^s1 y(t,r,s2=tau) dr;
% %             y_{t} = -y_{s2};
% % With BCs    v(t,s1=0) = 0;   e(t,s1=pi/2) = 0;   y(t,s1,s2=0) = e(t,s1).
% %
% % Parameter tau can be set; default tau = 0.9*tau*.
% %
% % Threshold [derived from the modes]: stable iff tau < tau* = 2.881401,
% % the delay margin of mode k = 1.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman
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

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2 r1

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
L = pi/2;
e = pde_var(s1,[0,L],1);
v = pde_var(s1,[0,L],1);
y = pde_var([s1;s2],[0,L;0,tau],[0;1]);         % no derivative in s1, 1st in s2
PDE_t = [diff(e,'t')==diff(v,s1);
         diff(v,'t')==0.1*diff(e,s1)-2*v-int(subs(subs(y,s2,tau),s1,r1),r1,[0,s1]);
         diff(y,'t')==-diff(y,s2);
         subs(v,s1,0)==0;           subs(e,s1,L)==0;
         subs(y,s2,0)==e];


% Example-specific settings (second output only): the light tier, with
if nargout>1
    st = lpisettings('light');
    st.lpi.ineq.dD = 0;             % slack lift degree at its minimum Dmin
    st.sop.deg_Q{2} = [1 2 2];      % Q of degree 1 (multiplier), 2 (integrals)
    st.lpi.ineq.psatz{2} = 'none';  % 2-D slack with the plain term only
    st.lpi.ineq.cap_inc = 'auto';   % slack integral bases: theta+s <= Dmax-1
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
