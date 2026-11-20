function [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Wave_Eq_Damped(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_2D_WAVE_EQ_DAMPED(GUI,PARAMS) returns the
% damped 2-D wave equation below, in the state u = [phi; phi_t], as a
% 'pde_struct' object, and recommended settings ST (the light tier). It
% sets no executive flag; a nonzero GUI prints that no GUI representation
% exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'kap=2;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         recommended settings (lpisettings_sop structure); the file
%               sets no executive flag.
%
% %---------------------------------------------------------------------% %
% % % 2D wave equation with stabilizing feedback (Jagt & Peet, 2025, Sec. 7.2.2,
% % % Eqns. (18)-(19)):
% % % PDE         phi_{tt} = phi_{s1s1} + phi_{s2s2} - 2*kap*phi_{t} - kap^2*phi,
% % % With BCs    phi(s1=0) = 0;   phi_{s1}(s1=1) = 0;
% % %             phi(s2=0) = 0;   phi_{s2}(s2=1) = 0,       (s1,s2) in [0,1]^2.
% % %
% % % Use the state u = [u1; u2] = [phi; phi_{t}], both in S_2^{(2,2)}, then
% % %
% % % PDE         u_{t} = [0, 1; -kap^2, -2*kap]*u + [0, 0; 1, 0]*(u_{s1s1} + u_{s2s2})
% % % With BCs    u(s1=0) = 0;   u_{s1}(s1=1) = 0;
% % %             u(s2=0) = 0;   u_{s2}(s2=1) = 0,
% % %
% % % i.e. u2 = phi_t satisfies the same BCs and regularity as u1 = phi, as in
% % % (19) of the paper.
% %
% % Parameter kap (>= 0) can be set. For kap > 0 the system is exponentially
% % PIE-to-PDE stable with maximal rate kap (paper, Appx. C.2): phi =
% % exp(-kap*t)*psi with psi an undamped wave; kap = 0 is the undamped wave
% % itself, which is not exponentially stable. The paper certifies the rate
% % with the LPI of its Cor. 35.
% %
% % Reference values (paper, Table 2): the maximal certified rate for
% % kap = 1..7 is 0.9999, 1.9998, 2.9874, 3.9973, 4.9980, 5.9906, 6.9645.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_2D_Wave_Eq_Damped
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
% Initial coding MMP - 09/28/2026

% Initialize variables
pvar s1 s2

% Default parameter, overridden by the assignment strings in params
kap = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% % % Construct the PDE.
clear stateNameGenerator
u = pde_var(2,[s1;s2],[0,1;0,1],[2;2]);   % u = [phi; phi_t], order (2,2) in both components
PDE_t = [diff(u,'t')==[0,1;-kap^2,-2*kap]*u + [0,0;1,0]*(diff(u,s1,2)+diff(u,s2,2));
         subs(u,s1,0)==0;     subs(diff(u,s1),s1,1)==0;
         subs(u,s2,0)==0;     subs(diff(u,s2),s2,1)==0];

% Example-specific settings: the light tier
st = lpisettings_sop('light');

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% D. S. Jagt and M. M. Peet, "A State-Space Representation of Coupled Linear
% Multivariate PDEs and Stability Analysis using SDP," arXiv:2508.14840, 2025.
