function [PDE_t, st] = PIETOOLS_PDE_Ex_Transport_Edge_Diffusion(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_TRANSPORT_EDGE_DIFFUSION(GUI,PARAMS) returns
% the transport equation on an edge driving a diffusion in the normal
% direction below, as a 'pde_struct' object, and (second output) the
% recommended settings ST of its stability test. It sets stability = 1 and
% stability_form = 'Q' (generalized storage T'Q) in the base workspace; a
% nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'r2=5;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         (optional) recommended settings (lpisettings
%               structure) for the stability test of this example.
%
% %---------------------------------------------------------------------% %
% % Transport on an edge driving a diffusion in the normal direction:
% % PDE         x1_{t}(t,s1)    = -x1_{s1};               s1 in [0,1]
% %             x2_{t}(t,s1,s2) = x2_{s2s2} + r2*x2;      (s1,s2) in [0,1]^2
% % With BCs    x1(t,s1=0) = 0;
% %             x2(t,s1,s2=0) = x1(t,s1);   x2(t,s1,s2=1) = 0.
% % x2 carries no derivative in s1 (order 0 in s1, 2 in s2).
% %
% % Parameter r2 can be set; default r2 = 0.9*r2*.
% %
% % Threshold [derived]: x1(t,.) = 0 for t >= 1 (transport with zero
% % inflow); from then on each slice x2(t,s1,.) obeys the Dirichlet heat
% % equation on [0,1] with reaction r2, with eigenvalues r2 - k^2*pi^2,
% % k >= 1:
% %   exponentially stable iff r2 < r2* = pi^2.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Transport_Edge_Diffusion
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

% Default parameter (0.9*r2*), overridden by the assignment strings in params
r2 = 0.9*pi^2;
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
x1 = pde_var(1,s1,[0,1]);
x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);       % order 0 in s1, 2 in s2
PDE_t = [diff(x1,'t')==-diff(x1,s1);
         diff(x2,'t')==diff(x2,s2,2)+r2*x2;
         subs(x1,s1,0)==0;
         subs(x2,s2,0)==x1;         subs(x2,s2,1)==0];


% Example-specific settings (second output only): the light tier, with
if nargout>1
    st = lpisettings('light');
    st.lpi.ineq.dD = 0;             % slack lift degree at its minimum Dmin
    st.sop.deg_Q{2} = 2;            % Q of degree 2
    st.lpi.ineq.wcross = 1;         % slack weight 1 in non-owned directions
end

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
