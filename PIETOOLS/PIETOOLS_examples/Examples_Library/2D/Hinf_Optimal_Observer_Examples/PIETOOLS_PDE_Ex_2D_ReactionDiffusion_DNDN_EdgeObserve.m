function [PDE_t,st] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_2D_REACTIONDIFFUSION_DNDN_EDGEOBSERVE(GUI,
% PARAMS) returns the 2-D reaction-diffusion equation below, with
% disturbance w1, sensed outputs y1, y2 with noise w2, w3, and integral
% output z, as a 'pde_struct' object, and the recommended settings ST of
% its H-infinity optimal estimator synthesis. It sets Hinf_estimator = 1
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
% - st:         example-specific settings (lpisettings structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Reaction-diffusion input-output-PDE with 
% %     Dirichlet-Neumann Dirichlet-Neumann BCs
% % PDE         x_{t}   = r*x + nu*x_{s1s1} + nu*x_{s2s2} + Cw*w1;
% % OUT         y1(s1)  = int_{c}^{d} x(t,s1,s2) ds2 + w2;
% %             y2(s2)  = int_{a}^{b} x(t,s1,s2) ds1 + w3;
% %             z       = int_{a}^{b} int_{c}^{d} x(t,s1,s2) ds2 ds1
% % With BCs    x(s1=a) = 0;        x(s2=c) = 0;
% %             x_{s1}(s1=b) = 0;   x_{s2}(s2=d) = 0;
% % (The commented term-based block below senses the edges instead:
% %  y1 = x(s2=d) + w2, y2 = x(s1=b) + w3.)
% %
% % Parameters nu and r, the state size ne, and limits a, b, c, d can be set.
% % Function Cw can also be set, defaults to
% %     Cw(s1,s2) = (s1^2-b)*(s2^2-d);
% % Stable for r <= mu, where
% %     mu = nu*pi^2*((1/4)/(b-a)^2 +(1/4)/(d-c)^2)
% %
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve
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
% Initial coding DJ - 08/02/2024
% MMP, 11/21/2026: Changes since PIETOOLS 2025: second output st, the
%   recommended settings of the H-infinity optimal estimator synthesis;
%   help revised, its outputs y1, y2 corrected to the integrals the code
%   takes; license added.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pde_struct PDE_t;
pvar s1 s2

%%% Executive Function:
evalin('base','Hinf_estimator = 1;');

% Default parameters, overridden by the assignment strings in params
nu = 1;     r = 8;   ne = 1;
a = 0;  b=1;        c = 0;       d = 1;
Cw = (s1^2-b)*(s2^2-d);
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
x = pde_var(ne,[s1;s2],[a,b;c,d]);
w1 = pde_var('in',ne);
w2 = pde_var('in',ne);
w3 = pde_var('in',ne);
y1 = pde_var('sense',ne,s1,[a,b]);
y2 = pde_var('sense',ne,s2,[c,d]);
z = pde_var('out',ne);
PDE_t = [diff(x,'t')==r*x+nu*(diff(x,s1,2)+diff(x,s2,2))+Cw*w1;
         y1==int(x,s2,[c,d])+w2;
         y2==int(x,s1,[a,b])+w3;
         z==int(x,[s1;s2],[a,b;c,d]);
         subs(x,s1,a)==0;       subs(diff(x,s1),s1,b)==0;
         subs(x,s2,c)==0;       subs(diff(x,s2),s2,d)==0];


% %%% Term-based input format
% PDE_t.x{1}.vars = [s1;s2];      PDE_t.x{1}.dom = [a,b;c,d];
% PDE_t.x{1}.size = ne;
% PDE_t.w{1}.vars = [];           PDE_t.w{1}.size = ne;
% PDE_t.w{2}.vars = [];           PDE_t.w{2}.size = ne;
% PDE_t.w{3}.vars = [];           PDE_t.w{3}.size = ne;
% PDE_t.y{1}.vars = s1;           PDE_t.y{1}.dom = [a,b];
% PDE_t.y{1}.size = ne;
% PDE_t.y{2}.vars = s2;           PDE_t.y{2}.dom = [c,d];
% PDE_t.y{2}.size = ne;
% PDE_t.z{1}.vars = [];           PDE_t.z{1}.size = ne;
% 
% % PDE: x_{t} = [r, c1, c2] * [x; x_{s1s1}; x_{s2s2}] + (s1^2-1)*(s2^2-1)*w
% PDE_t.x{1}.term{1}.D = [0,0; 2,0; 0,2];
% PDE_t.x{1}.term{1}.C = [r*eye(ne), nu*eye(ne), nu*eye(ne)];
% 
% PDE_t.x{1}.term{2}.w = 1;
% PDE_t.x{1}.term{2}.C = Cw;
% 
% % OUT: y2(t,s2) = x(t,s1=b,s2) + w3;
% PDE_t.y{2}.term{1}.x = 1;   PDE_t.y{2}.term{1}.loc = [b,s2];
% PDE_t.y{2}.term{2}.w = 3;
% 
% % OUT: y1(t,s1) = x(t,s1,s2=d) + w2;
% PDE_t.y{1}.term{1}.x = 1;   PDE_t.y{1}.term{1}.loc = [s1,d];
% PDE_t.y{1}.term{2}.w = 2;
% 
% % OUT: z(t) = int_{a}^{b} int_{c}^{d} x(t,s1,s2) ds2 ds1
% PDE_t.z{1}.term{1}.x = 1;
% PDE_t.z{1}.term{1}.I = {[a,b];[c,d]};
% 
% % BC1: 0 = x(s1,c)                     % BC3: 0 = x_{s2}(s1,d)
% PDE_t.BC{1}.term{1}.loc = [s1,c];      PDE_t.BC{3}.term{1}.loc = [s1,d];
%                                        PDE_t.BC{3}.term{1}.D = [0,1];
% % BC2: 0 = x(a,s2)                     % BC4: 0 = x_{s1}(b,s2)
% PDE_t.BC{2}.term{1}.loc = [a,s2];      PDE_t.BC{4}.term{1}.loc = [b,s2];
%                                        PDE_t.BC{4}.term{1}.D = [1,0];


% Example-specific settings: the heavy tier                                 % MMP, 11/21/2026
st = lpisettings('heavy');                                                  % MMP, 11/21/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end