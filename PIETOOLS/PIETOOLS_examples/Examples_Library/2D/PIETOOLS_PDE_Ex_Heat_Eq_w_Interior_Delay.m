function [PDE_t, st] = PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay(GUI,params) % MMP, 11/17/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_HEAT_EQ_W_INTERIOR_DELAY(GUI,PARAMS) returns
% the reaction-diffusion equation below with a delayed reaction term, the
% delay written as a transport state x2(s1,s2), as a 'pde_struct' object,
% and (second output) the recommended settings ST of its stability test.
% It sets stability = 1 and stability_form = 'Q' (generalized storage T'Q)
% in the base workspace; a nonzero GUI prints that no GUI representation
% exists.
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
% - st:         (optional) recommended settings (lpisettings_sop            % MMP, 11/17/2026
%               structure) for the stability test of this example,          % MMP, 11/17/2026
%               returned only when requested (nargout > 1).                 % MMP, 11/17/2026
%
% %---------------------------------------------------------------------% %
% % Diffusive Equation with Delay in Dynamics (Caliskan 2009):
% % PDE         x_{t}  = c*x_{s1s1}(t,s1) + a0*x(t,s1) - a1*x(t-tau,s1);    s1 in [0,pi]   
% % With BCs    x(t,s1=0) = 0;   
% %             x(t,s1=pi) = 0;
% %
% % Let x1=x, and x2(t,s1,s2)=x1(t-s2,s1), introducing transport equation:
% %         x1(t) = x2(t,0),   x2_{t} = -x2_{s2}        s2 in [0,tau]
% % Then
% %
% % PDE:    x1_{t}         = c*x1_{s1s1}(t,s1) + a0*x1(t,s1) - a1*x2(t,s1,s2=tau);
% %         x2_{t}         = -x2_{s2};
% % BCs:    x1(t,s1=0)     = 0;
% %         x1(t,s1=pi)    = 0;
% %         x2(t,s1,s2=0)  = x1(t,s1);
% %         x2(t,s1=0,s2)  = 0;
% %         x2(t,s1=pi,s2) = 0;
% %
% % Parameters tau, c, a0, a1 and the domain [a,b] (default [0,pi]) can be set.
% % For c = 1, a0 = 1.9, and a1 = 1, stable for tau<1.0347.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay
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
% Initial coding DJ - 11/30/2022
% MMP, 11/17/2026: Second output st, the settings that certify the
%                   example (maintainer: application-specific settings
%                   where the tiers do not); the light tier certifies it
%                   in neither form (header);
% function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay(GUI,params)   % MMP, 11/17/2026 (was)
% MMP, 11/18/2026: st: the slack lift degree at its minimum and Q of
%                   degree 2 in two variables (deg_Q), in place of the
%                   dw_space raise of 11/17: half the variables (header);
% MMP, 11/19/2026: st updated to the per-example search of 11/19/2026
%                   (maintainer: no default settings): see the block;
%                   the form of its tests named in the workspace.
% MMP, 11/19/2026: the settings notes of the header (11/17-11/19) removed: they
%                   recorded tests and files the user does not have, not the problem;
% MMP, 11/19/2026: the remaining settings notes removed: they recorded tests,
%                   not the problem;
% MMP, 11/19/2026: the dated test records of the header removed: they
%                   recorded tests and files the user does not have, not
%                   the problem;
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026
%evalin('base','stability_dual = 1;');

% Default parameters, overridden by the assignment strings in params
c = 1;      a0 = 1.9;       a1 = 1;
ne = 1;     tau = 1; 
a = 0;      b = pi;
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
x1 = pde_var(s1,[a,b]);
x2 = pde_var([s1;s2],[a,b;0,tau],[2;1]);    % x2 must be 2nd order differentiable wrt s1, and 1st order wrt s2
PDE_t = [diff(x1,'t')==c*diff(x1,s1,2)+a0*x1-a1*subs(x2,s2,tau);
         diff(x2,'t')==-diff(x2,s2);
         subs([x1;x2],s1,a)==0;     subs([x1;x2],s1,b)==0;
         subs(x2,s2,0)==x1];


% %%% Term-based input format
% % Initialize the state variables and input.
% pde_struct PDE_t;
% PDE_t.x{1}.vars = s1;           PDE_t.x{1}.dom = [a,b];
% PDE_t.x{2}.vars = [s1;s2];      PDE_t.x{2}.dom = [a,b;0,tau];
% PDE_t.x{2}.diff = [2,1];    % x2 must be 2nd order differentiable wrt s1, and 1st order wrt s2
% 
% % PDE: x1_{t}(t,s1) = c*x1_{s1s1}(t,s1)
% PDE_t.x{1}.term{1}.x = 1;
% PDE_t.x{1}.term{1}.D = 2;
% PDE_t.x{1}.term{1}.C = c;
% 
% % PDE: x1_{t}(t,s1) = ... + a0*x1(t,s1)
% PDE_t.x{1}.term{2}.x = 1;
% PDE_t.x{1}.term{2}.C = a0;
% 
% % PDE: x1_{t}(t,s1) = ... + a1*x2(t,s1,tau)
% PDE_t.x{1}.term{3}.x = 2;
% PDE_t.x{1}.term{3}.loc = [s1,tau];
% PDE_t.x{1}.term{3}.C = -a1;
% 
% % PDE: x2_{t}(t,s1,s2) = -x2_{s2}(t,s1,s2)
% PDE_t.x{2}.term{1}.x = 2;
% PDE_t.x{2}.term{1}.D = [0,1];
% PDE_t.x{2}.term{1}.C = -1;
% 
% 
% % BC1: x1(s1=a) = 0;                BC2: x1(s1=b)=0;
% PDE_t.BC{1}.term{1}.x = 1;          PDE_t.BC{2}.term{1}.x = 1;
% PDE_t.BC{1}.term{1}.loc = a;        PDE_t.BC{2}.term{1}.loc = b;
% 
% % BC3: (Sewing condition) x2(s1,s2=0) = x1(s1)
% PDE_t.BC{3}.term{1}.x = 2;          PDE_t.BC{3}.term{2}.x = 1;
% PDE_t.BC{3}.term{1}.loc = [s1,0];   PDE_t.BC{3}.term{2}.C = -1;
% 
% % BC4: x2(s1=a,s2) = 0;             BC5: x2(s1=b,s2) = 0;
% PDE_t.BC{4}.term{1}.x = 2;          PDE_t.BC{5}.term{1}.x = 2;
% PDE_t.BC{4}.term{1}.loc = [a,s2];   PDE_t.BC{5}.term{1}.loc = [b,s2];


% Application-specific settings (header): the light tier, the slack         % MMP, 11/18/2026 (was); MMP, 11/17/2026
% weight of the space of x1 (L2[s1], the first space) raised by one         % MMP, 11/18/2026 (was); MMP, 11/17/2026
% Example-specific settings                                                 % MMP, 11/19/2026
% (second output only; light tier; slack lift degree at its minimum Dmin;
% Q of degree 1 in the multiplier, 2 in the integrals; slack integral
% bases capped at total degree theta+s <= Dmax-1)
if nargout>1                                                                % MMP, 11/17/2026
    st = lpisettings_sop('light');                                          % MMP, 11/17/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/18/2026
%   st.sop.deg_Q{2} = 2;                                                    % MMP, 11/19/2026 (was); MMP, 11/18/2026
    st.sop.deg_Q{2} = [1 2 2];                                              % MMP, 11/19/2026
%    st.lpi.ineq.dw_space = [1 0];                                          % MMP, 11/18/2026 (was); MMP, 11/17/2026
    st.lpi.ineq.cap_inc = 'auto';                                           % MMP, 11/19/2026
end                                                                         % MMP, 11/17/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{CALISKAN2009,
% title = {Stability Analysis of the Heat Equation with Time-Delayed Feedback},
% journal = {IFAC Proceedings Volumes},
% volume = {42},
% number = {6},
% pages = {220-224},
% year = {2009},
% note = {6th IFAC Symposium on Robust Control Design},
% issn = {1474-6670},
% doi = {https://doi.org/10.3182/20090616-3-IL-2002.00038},
% url = {https://www.sciencedirect.com/science/article/pii/S1474667015404057},
% author = {Sina Yamaç çalişkan and Hitay özbay},
% }