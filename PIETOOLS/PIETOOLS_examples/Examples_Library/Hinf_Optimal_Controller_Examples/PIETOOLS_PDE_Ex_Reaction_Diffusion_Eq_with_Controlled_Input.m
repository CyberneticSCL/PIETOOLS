function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Controlled_Input(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_REACTION_DIFFUSION_EQ_WITH_CONTROLLED_
% INPUT(GUI,PARAMS) returns the reaction-diffusion equation below, with
% distributed control u, disturbance w and regulated output z, in the
% term-based and in the batch format, and the recommended settings ST of
% its H-infinity optimal controller synthesis. It sets Hinf_control = 1 in
% the base workspace. With GUI = 1 it also loads the associated GUI save
% file into a new PIETOOLS_PDE_GUI window.
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Reaction-diffusion equation with distributed control u:
% % PDE                            x_{t} = lam*x + x_{ss} + (s-s^2)u(t) + (s-s^2)w
% % With BCs                         x(s=0) = 0
% %                                         x(s=1) = 0
% % and regulated output  z =[ int(x(s,t),s,0,1);u]
% % Parameter lam can be set.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Controlled_Input
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
% Initial coding DJ - 08/18/2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: third output st, the
%   recommended settings of the H-infinity optimal controller synthesis;
%   help revised (title of the system corrected), license added.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Controlled_Input(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_control = 1;')

% Default parameter, overridden by the assignment strings in params
lam = 10;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Batch input format
% number of inputs
PDE_b.nw = 1;   PDE_b.nu = 1;
% number of outputs
PDE_b.ny = 0;   PDE_b.nz = 2; 
% number of ODE states
PDE_b.nx = 0;
% number of PDE states
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;
PDE_b.dom = [0,1];
% PDE state dynamics
PDE_b.A0 = lam;   PDE_b.A2 = 1;
% control input to PDE state
PDE_b.B22 = s-s^2;   
% disturbance to PDE state
PDE_b.B21 = s-s^2;    
% boundary conditions
PDE_b.B = [1 0 0 0;
           0 1 0 0];
% PDE state to regulated output
PDE_b.Ca1 = [ 1
                        0];
% control input to regulated output
PDE_b.D12 = [0
                        1];


%%% Term-based input format
% Initialize 1D PDE state component.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
% Initialize finite-dimensional inputs.
PDE_t.u{1}.vars = [];
PDE_t.w{1}.vars = [];

% Initialize finite-dimensional, vector-valued outputs
PDE_t.z{1}.vars=[];
PDE_t.z{2}.vars=[];

% PDE: x_{t} = lam * x + x_{ss}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = [0; 2];
PDE_t.x{1}.term{1}.C = [lam, 1];

% PDE: x_{t} = ... + w
PDE_t.x{1}.term{2}.w = 1;
PDE_t.x{1}.term{2}.C = s-s^2;
% PDE: x_{t} = ... + u
PDE_t.x{1}.term{3}.u = 1;
PDE_t.x{1}.term{3}.C = s-s^2;
% BC 1: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x(1)
PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.loc = 1;

% z = [z_1 ; z_2] such that 
% z_1 = \int_0^1 x(s,t) ds
PDE_t.z{1}.term{1}.x = 1;
PDE_t.z{1}.term{1}.I{1} = [0, 1];
% z_2 = u
PDE_t.z{2}.term{1}.u = 1;
% command line
% pvar t s;
% x = state('pde'); u = state('in');
% w = state('in'); z=state('out',2);
% pde = sys();
% eq_dyn = diff(x,t)==10*x+diff(x,s,2)+(s-s^2)u+(s-s^2)w;
% eq_out = z==[int(x,s,[0,1]);u];
% pde = addequation(pde,[eq_dyn;eq_out]);
% eq_bc = [subs(x, s, 0) == 0;subs(x, s, 1) == 0];
% pde = addequation(pde,eq_bc);

if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'Examples_Libary/Hinf_Optimal_Controller_Examples/Reaction_Diffusion_Eq_with_Controlled_Input.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/21/2026
% (stripped tier; 1-D storage from the lower-integral basis operator only,
% theta-length 2: the multiplier and the upper integral have length 0;
% slack integral bases capped at total degree theta+s <= Dmax-1; free
% operators of degree 0 in the multiplier and 2 in the integral kernels)
st = lpisettings_sop('stripped');                                           % MMP, 11/21/2026
st.lpi.pos.len_com{1} = 0;                                                  % MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [0;2;2];                                            % MMP, 11/21/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/21/2026
st.lpi.var.deg{1} = [0 2 2];                                                % MMP, 11/21/2026

end