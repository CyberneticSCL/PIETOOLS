function [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_HEAT_DELAYED_ROBIN_BOUNDARY_NICAISE(GUI,
% PARAMS) returns the heat equation with a delayed Robin boundary term
% below, the delayed boundary value written as a transport state y(s2), as
% a 'pde_struct' object, and the recommended settings ST of its stability
% test. It sets stability = 1 and stability_form = 'Q' (generalized
% storage T'Q) in the base workspace; a nonzero GUI prints that no GUI
% representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'mu1=2;','tau=1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         recommended settings (lpisettings structure) for the
%               test this file names.
%
% %---------------------------------------------------------------------% %
% % Heat equation with a delayed Robin boundary term (Nicaise, Valein,
% % Fridman 2009, system (1), constant delay, a = 1):
% % PDE         u_{t}(t,s1) = u_{s1s1}(t,s1);                s1 in [0,pi]
% % With BCs    u(t,s1=0) = 0;
% %             u_{s1}(t,pi) = -mu0*u(t,pi) - mu1*u(t-tau,pi);
% %
% % Let y(t,s2) = u(t-tau*s2,pi), s2 in [0,1] (system (6) of the paper):
% % PDE         u_{t} = u_{s1s1};   y_{t} = -(1/tau)*y_{s2};
% % With BCs    u(t,s1=0) = 0;
% %             u_{s1}(t,s1=pi) = -mu0*u(t,pi) - mu1*y(t,s2=1);
% %             y(t,s2=0) = u(t,s1=pi);
% %
% % Parameters mu0, mu1 and tau can be set; default mu0 = 1, mu1 = 0.5,
% % tau = 1 (values chosen here; the paper gives none).
% %
% % The paper gives sufficient conditions for exponential stability only
% % [abstract]; no threshold is known here. Spectral abscissa by the method
% % of lines (200 heat and 400 transport nodes) [finite differences]:
% %   -0.7623 for mu1 = 0.5;   -0.1987 for mu1 = 2.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise
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
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test

% Default parameters, overridden by the assignment strings in params
mu0 = 1;    mu1 = 0.5;  tau = 1;
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
u = pde_var(s1,[0,pi]);
y = pde_var(s2,[0,1],1);
PDE_t = [diff(u,'t')==diff(u,s1,2);
         diff(y,'t')==-(1/tau)*diff(y,s2);
         subs(u,s1,0)==0;
         subs(diff(u,s1),s1,pi)==-mu0*subs(u,s1,pi)-mu1*subs(y,s2,1);
         subs(y,s2,0)==subs(u,s1,pi)];


% Example-specific settings: the stripped tier, with the slack lift degree
% at its minimum Dmin, the 2-D slack with the plain term only (psatz
% 'none'), the slack integral bases capped at theta+s <= Dmax-1, and Q of
% degree 2 in the multiplier and 1 in the integral kernels
st = lpisettings('stripped');
st.lpi.ineq.dD = 0;
st.lpi.ineq.psatz{2} = 'none';
st.lpi.ineq.cap_inc = 'auto';
st.sop.deg_Q{2} = [2 1 1];

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{nicaise2009stability,
%   title={Stability of the heat and of the wave equations with boundary time-varying delays},
%   author={Nicaise, Serge and Valein, Julie and Fridman, Emilia},
%   journal={Discrete and Continuous Dynamical Systems - Series S},
%   volume={2},
%   number={3},
%   pages={559--581},
%   year={2009},
%   doi={10.3934/dcdss.2009.2.559}
% }
