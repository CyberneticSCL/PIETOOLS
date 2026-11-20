function [PDE_t,st] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDN(GUI,params) % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_2D_REACTIONDIFFUSION_DDDN(GUI,PARAMS)
% returns the 2-D reaction-diffusion equation below, with disturbance w,
% integral output z, Dirichlet BCs on three edges and a Neumann BC at
% s2 = d, as a 'pde_struct' object, and the recommended settings ST of its
% H-infinity gain test. It sets Hinf_gain = 1 and Hinf_gain_form = 'P'
% (storage <T x, R T x>, R >= eppos I) in the base workspace; a nonzero GUI
% prints that no GUI representation exists.
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Reaction-diffusion input-output-PDE with 
% %     Dirichlet-Dirichlet Dirichlet-Neumann BCs
% % PDE         x_{t}   = r*x + nu*x_{s1s1} + nu*x_{s2s2} + w
% % OUT         z       = int_{a}^{b} int_{c}^{d} x(t,s1,s2) ds2 ds1
% % With BCs    x(s1=a) = 0;   x(s2=c) = 0;
% %             x(s1=b) = 0;   x_{s2}(s2=d) = 0;
% %
% % Parameters nu and r, the state size ne, and limits a, b, c, d can be set.
% % Stable for r <= mu_{1,1}, where
% %     mu_{m,n} = nu*pi^2*(m^2/(b-a)^2 +(n-1/2)^2/(d-c)^2)
% % Then, the L2-gain given by
% % gam = 8*sqrt((b-a)*(d-c))/pi^2 
% %             *sqrt(sum_{m=1}^{infty}sum_{n=1}^{infty} 1/((mu_{2m-1,n}-r)*(2m-1)*(2n-1))^2
% %
% % For nu = 1, r = 5, [a,b]=[c,d]=[0,1], we have gam = 0.11101;
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDN
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
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDN(GUI,params)  % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pde_struct PDE_t;
pvar s1 s2

%%% Executive Function:
evalin('base','Hinf_gain = 1;');
evalin('base','Hinf_gain_form = ''P'';');  % form of the gain test          % MMP, 11/19/2026
%evalin('base','Hinf_gain_dual = 1;');

% Default parameters, overridden by the assignment strings in params
nu = 1;     r = 5;   ne = 1;
a = 0;  b = 1;       c = 0;  d = 1;
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
w = pde_var('in',ne,[s1;s2],[a,b;c,d]);
z = pde_var('out',ne);
PDE_t = [diff(x,'t')==r*x+nu*(diff(x,s1,2)+diff(x,s2,2))+w;
         z==int(x,[s1;s2],[a,b;c,d]);
         subs(x,s1,a)==0;       subs(x,s1,b)==0;
         subs(x,s2,c)==0;       subs(diff(x,s2),s2,d)==0];


% %%% Term-based input format
% PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [a,b;c,d];
% PDE_t.x{1}.size = ne;
% PDE_t.w{1}.vars = [s1;s2];   PDE_t.w{1}.dom = [a,b;c,d];
% PDE_t.w{1}.size = ne;
% PDE_t.z{1}.vars = [];
% PDE_t.z{1}.size = ne;
% 
% % PDE: x_{t} = [r, c1, c2] * [x; x_{s1s1}; x_{s2s2}] + w
% PDE_t.x{1}.term{1}.D = [0,0; 2,0; 0,2];
% PDE_t.x{1}.term{1}.C = [r*eye(ne), nu*eye(ne), nu*eye(ne)];
% 
% PDE_t.x{1}.term{2}.w = 1;
% 
% % OUT: z(t) = int_{a}^{b} int_{c}^{d} x(t,s1,s2) ds2 ds1
% PDE_t.z{1}.term{1}.x = 1;
% PDE_t.z{1}.term{1}.I = {[a,b];[c,d]};
% 
% % BC1: 0 = x(s1,c)                     % BC3: 0 = x_{s2}(s1,d)
% PDE_t.BC{1}.term{1}.loc = [s1,c];      PDE_t.BC{3}.term{1}.loc = [s1,d];
%                                        PDE_t.BC{3}.term{1}.D = [0,1];
% % BC2: 0 = x(a,s2)                     % BC4: 0 = x(b,s2)
% PDE_t.BC{2}.term{1}.loc = [a,s2];      PDE_t.BC{4}.term{1}.loc = [b,s2];


% Example-specific settings                                                 % MMP, 11/19/2026
% (heavy tier; storage on the 2-D space with theta-lengths [2 2] for the
% multiplier and [1 0] for the integral, so no integral in s2; slack
% weight raise of the gain LPI 1)
st = lpisettings_sop('heavy');                                              % MMP, 11/19/2026
st.lpi.pos.len_sos{2} = [2 2;1 0;0 0];                                      % MMP, 11/19/2026
st.sop.dw_P = 1;                                                            % MMP, 11/19/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end


%%
% For comparison, compute the L2-gain (approximately)
M = 10;     N = 10;
Lx = b-a;   Ly = d-c;
r_lim = nu*pi^2*(1/Lx^2 +(1/4)/Ly^2);       % Stable iff r<=r_lim
mu_mn = nu*pi^2*(((2*(1:M)'-1)/Lx).^2 +(((1:N)-(1/2))/Ly).^2);
mn_fact = ((2*(1:M)'-1)).*((2*(1:N)-1));
full_fact = 1./((mu_mn-r).*mn_fact);
fct = sqrt(sum(full_fact.^2,'all'));
L2_gain = (8*sqrt(Lx*Ly)/(pi^2))*fct;       % L2 gain (approximately)

end