function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Robin(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_REACTION_DIFFUSION_ROBIN(GUI,PARAMS)
% returns the reaction-diffusion equation below with Robin BCs, in the
% pde_var (term-based) and in the batch format, and the recommended
% settings ST of its stability tests. It sets stability = 1,
% stability_dual = 1 and stability_form = 'P' (storage <T x, P T x>,
% P >= eppos I; dual form 'Pd') in the base workspace; a nonzero GUI
% prints that no GUI representation exists. The comment after the
% function derives the stability threshold.
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
% - PDE_b:      PDE structure defining the example system in the batch
%               format.
% - st:         example-specific settings (lpisettings structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Reaction-Diffusion Equation with Robin Boundary Conditions:
% % PDE        x_{t} = lam*x + x_{ss} 
% % with BCs   x(s=0) = c1*x_{s}(s=0)
% %            x(s=1) = -c2*x_{s}(s=1)
% % 
% % Parameters lam, c1, c2, and ne (state size) can be set.
% % When c1=c2=1, unstable for lam > 1.7071.
% % Default parameters lam=1.705, c1=c2=1, ne=1.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Reaction_Diffusion_Robin
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
% MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
%   recommended settings of the stability tests; stability_form = 'P' set
%   in the base workspace; help and comments revised.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Robin(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;')
evalin('base','stability_form = ''P'';');  % form of the stability test     % MMP, 11/21/2026

% Specify the parameters
ne = 1;   lam = 1.705; %1.705
c1 = 1;     c2 = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end



% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2=ne;
PDE_b.dom = [0,1];

PDE_b.A0 = lam*eye(ne);
PDE_b.A2 = 1*eye(ne);

PDE_b.B = [eye(ne) -c1*eye(ne) zeros(ne) zeros(ne);
           zeros(ne)  zeros(ne) eye(ne) c2*eye(ne)];


%%% pde_var input format
clear stateNameGenerator
x = pde_var(ne,s,[0,1]);
PDE_t = [diff(x,'t')==lam*x+diff(x,s,2);
         subs(x,s,0)==c1*subs(diff(x,s),s,0);
         subs(x,s,1)==-c2*subs(diff(x,s),s,1)];


% %%% Term-based input format
% PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
% 
% % PDE: x_{t} = lam*x
% PDE_t.x{1}.term{1}.C = lam*eye(ne);
% 
% % PDE: x_{t} = ... + x_{ss}
% PDE_t.x{1}.term{2}.D = 2;
% 
% % BCs: 0 = x(s=0) - c1*x_{s}(s=0)
% PDE_t.BC{1}.term{1}.loc = 0;        PDE_t.BC{1}.term{2}.loc = 0;
%                                     PDE_t.BC{1}.term{2}.D = 1;
%                                     PDE_t.BC{1}.term{2}.C = -c1*eye(ne);
% 
% % BCs: 0 = x(s=1) + c2*x_{s}(s=1)
% PDE_t.BC{2}.term{1}.loc = 1;        PDE_t.BC{2}.term{2}.loc = 1;
%                                     PDE_t.BC{2}.term{2}.D = 1;
%                                     PDE_t.BC{2}.term{2}.C = c2*eye(ne);


% No GUI representation of this system
if GUI
    disp('No GUI representation available for this system.')
end

% Example-specific settings                                                 % MMP, 11/21/2026
% (stripped tier; 1-D storage with upper-integral s-length 1 and
% theta-lengths 0 (multiplier: absent), 1 (upper), 1 (lower integral))
st = lpisettings('stripped');                                               % MMP, 11/21/2026
st.lpi.pos.len_com{1} = 1;                                                  % MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [0;1;1];                                            % MMP, 11/21/2026

end

%% Analytic stability result:
% We solve the equation, using separation of variables, x(t,s)=T(t)S(s)
% -->   T'(t)*S(s) = lam*T(t)S(s)+T(t)S''(s)
% -->   T'(t)/T(t) = lam + S''(s)/S(s)
% -->   T'(t)/T(t) = lam-mu     S''(s)/S(s) = -mu
% We get
%   T'(t) = lam-mu*T(t), so T(t) = exp((lam-mu)*t)*T(0)
% and
%   S''(s) = -mu*S(s), so for mu>0,
%       S(s) = a*sin(sqrt(mu)*s) + b*cos(sqrt(mu)*s)
% Then S(0) = b and S'(0)=a*sqrt(mu) so we must have b=c1*a*sqrt(mu)
% Furthermore, S(1)=a*sin(sqrt(mu))+b*cos(sqrt(mu))
%                S'(1)=a*sqrt(mu)*cos(sqrt(mu)) - b*sqrt(mu)*sin(sqrt(mu))
% so we must have 
% 0 = S(1)+c2*S'(1) 
%   = a*[sin(sqrt(mu))+c2*sqrt(mu)*cos(sqrt(mu))]
%           + b*[cos(sqrt(mu)) -c2*sqrt(mu)*sin(sqrt(mu)]
%   = a*[sin(sqrt(mu))+c2*sqrt(mu)*cos(sqrt(mu))]
%           + c1*a*sqrt(mu)*[cos(sqrt(mu)) -c2*sqrt(mu)*sin(sqrt(mu)]
%   = a*[(1-c1*c2*mu)*sin(sqrt(mu))+(c1+c2)*sqrt(mu)*cos(sqrt(mu))]
% Thus, we must have
%   (1-c1*c2*mu)*sin(sqrt(mu)) = -(c1+c2)*sqrt(mu)*cos(sqrt(mu))
% For c1=c2=1, we can numerically solve this equation, to find the smallest
% value sqrt(mu) for which this is satisfied to be sqrt(mu)=1.306542374,
% yielding mu = 1.7071. Since T(t) = exp((lam-mu)*t)*T(0), stability
% thus requires lam <= 1.7071, when c1=c2=1.
