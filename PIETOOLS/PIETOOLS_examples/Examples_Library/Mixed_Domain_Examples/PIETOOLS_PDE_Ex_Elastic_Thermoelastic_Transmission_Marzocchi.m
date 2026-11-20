function [PDE_t,st] = PIETOOLS_PDE_Ex_Elastic_Thermoelastic_Transmission_Marzocchi(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_ELASTIC_THERMOELASTIC_TRANSMISSION_MARZOCCHI(
% GUI,PARAMS) returns the 1-D elastic-thermoelastic transmission problem
% below, both parts written on s1 in [0,L], as a 'pde_struct' object, and
% the recommended settings ST of its stability test. It sets stability = 1
% and stability_form = 'P' (storage <T x, P T x>, P >= eppos I) in the
% base workspace; a nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'L=0.5;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         recommended settings (lpisettings_sop structure) for the
%               test this file names.
%
% %---------------------------------------------------------------------% %
% % Elastic-thermoelastic transmission problem, 1D (Marzocchi, Munoz
% % Rivera, Naso 2002) [unverified standard form], a = b = k = alpha =
% % beta = 1:
% % PDE         u_{tt} = u_{xx}                         x in [0,L]
% %             v_{tt} = v_{xx} - theta_{x},
% %             theta_{t} = theta_{xx} - v_{xt}         x in [L,2L]
% % Transmission at x = L: continuity of displacement and stress.
% %
% % Let e1 = u_{x}, p = u_{t} on [0,L], and e2 = v_{x}, q = v_{t}, th =
% % theta on [L,2L], shifted to s1 in [0,L]:
% % PDE         e1_{t} = p_{s1};   p_{t} = e1_{s1};
% %             e2_{t} = q_{s1};   q_{t} = e2_{s1} - th_{s1};
% %             th_{t} = th_{s1s1} - q_{s1};
% % With BCs    e1(t,s1=0) = 0;   q(t,s1=L) = 0;   th(t,0) = th(t,L) = 0;
% %             p(t,s1=L) = q(t,s1=0);   e1(t,s1=L) = e2(t,s1=0);
% % Modification: a free end at x = 0 (e1(0) = 0) replaces the clamped end
% % of the paper, which leaves the uniform-strain mode as an equilibrium of
% % the Riemann form (a zero eigenvalue of the PIE).
% %
% % Parameter L can be set; default L = 0.5.
% %
% % The paper shows exponential decay of the energy for the clamped
% % problem [secondary citations and abstract]. Exponential stability of
% % the modified system is not established here; no threshold.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Elastic_Thermoelastic_Transmission_Marzocchi
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
pvar s1

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test

% Default parameter, overridden by the assignment strings in params
L = 0.5;
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
e1 = pde_var(s1,[0,L],1);   pp = pde_var(s1,[0,L],1);
e2 = pde_var(s1,[0,L],1);   q = pde_var(s1,[0,L],1);
th = pde_var(s1,[0,L]);
PDE_t = [diff(e1,'t')==diff(pp,s1);
         diff(pp,'t')==diff(e1,s1);
         diff(e2,'t')==diff(q,s1);
         diff(q,'t')==diff(e2,s1)-diff(th,s1);
         diff(th,'t')==diff(th,s1,2)-diff(q,s1);
         subs(e1,s1,0)==0;          subs(q,s1,L)==0;
         subs(th,s1,0)==0;          subs(th,s1,L)==0;
         subs(pp,s1,L)==subs(q,s1,0);
         subs(e1,s1,L)==subs(e2,s1,0)];


% Example-specific settings: the stripped tier, with the 1-D storage from
% the lower-integral basis operator only (multiplier theta-length 0,
% upper-integral s-length 0), theta-length 1; slack lift degree at its
% minimum Dmin
st = lpisettings_sop('stripped');
st.lpi.pos.len_com{1} = 0;
st.lpi.pos.len_sos{1} = [0;2;1];
st.lpi.ineq.dD = 0;

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{marzocchi2002asymptotic,
%   title={Asymptotic behaviour and exponential stability for a transmission problem in thermoelasticity},
%   author={Marzocchi, Alfredo and Mu{\~n}oz Rivera, Jaime E and Naso, Maria Grazia},
%   journal={Mathematical Methods in the Applied Sciences},
%   volume={25},
%   number={11},
%   pages={955--980},
%   year={2002},
%   doi={10.1002/mma.323}
% }
