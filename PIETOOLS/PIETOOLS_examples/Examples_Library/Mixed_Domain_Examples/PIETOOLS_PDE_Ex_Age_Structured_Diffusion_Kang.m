function [PDE_t, st] = PIETOOLS_PDE_Ex_Age_Structured_Diffusion_Kang(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_AGE_STRUCTURED_DIFFUSION_KANG(GUI,PARAMS)
% returns the age-structured diffusion equation below, with a renewal
% (nonlocal) BC in age, as a 'pde_struct' object, and the recommended
% settings ST of its stability test. It sets stability = 1 and
% stability_form = 'Q' (generalized storage T'Q) in the base workspace; a
% nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'b=1;','d=0.1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         (optional) recommended settings (lpisettings_sop
%               structure) for the stability test of this example.
%
% %---------------------------------------------------------------------% %
% % Age-structured population with spatial diffusion (linear McKendrick
% % model with diffusion; Kang 2022) [standard linear form, unverified],
% % s1 = age a, s2 = space x:
% % PDE         q_{t}(t,s1,s2) = -q_{s1} + d*q_{s2s2} - mu*q;
% %                                               (s1,s2) in [0,1]^2
% % With BCs    q(t,s1=0,s2) = b*int_0^1 q(t,s1,s2) ds1;   (renewal)
% %             q(t,s1,s2=0) = 0;   q(t,s1,s2=1) = 0;
% % q is first order in s1 and second order in s2.
% %
% % Parameters b, d and mu can be set; default d = 0.1, mu = 1,
% % b = 0.9*b*.
% %
% % Threshold [derived from the modes, exact]: in the sine modes in s2 the
% % first mode decays iff R_eff = b*(1-exp(-c))/c < 1, c = mu + d*pi^2:
% %   stable iff b < b* = c/(1-exp(-c)) = 2.302685   (c = 1.986960).
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Age_Structured_Diffusion_Kang
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

% Default parameters (b = 0.9*b*), overridden by the strings in params
d = 0.1;    mu = 1;     b = 0.9*2.302685;
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
q = pde_var([s1;s2],[0,1;0,1],[1;2]);           % 1st order in s1 (age), 2nd in s2
PDE_t = [diff(q,'t')==-diff(q,s1)+d*diff(q,s2,2)-mu*q;
         subs(q,s1,0)==b*int(q,s1,[0,1]);
         subs(q,s2,0)==0;           subs(q,s2,1)==0];


% Example-specific settings (second output only): the light tier, with
if nargout>1
    st = lpisettings_sop('light');
    st.lpi.ineq.dD = 0;             % slack lift degree at its minimum Dmin
    st.sop.deg_Q{2} = [0 1 1];      % Q of degree 0 (multiplier), 1 (integrals)
end

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{kang2022effects,
%   title={Effects of diffusion on the principal eigenvalue for age-structured models with random diffusion},
%   author={Kang, Hao},
%   journal={Proceedings of the Royal Society of Edinburgh Section A},
%   volume={152},
%   number={1},
%   pages={258--280},
%   year={2022},
%   doi={10.1017/prm.2021.10}
% }
