function [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_HEAT_EQ_W_DELAYED_BOUNDARY_INPUT(GUI,PARAMS)
% returns the reaction-diffusion equation below with a delayed control
% input u(t-tau) at s1 = 1, as a 'pde_struct' object (the delay kept as a
% delay; the transport form below is the alternative branch of the code),
% and the recommended settings ST of its stability test. It sets
% stability = 1 and stability_form = 'Q' (generalized storage T'Q) in the
% base workspace; a nonzero GUI prints that no GUI representation exists.
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
% % Diffusive Equation with Delayed Boundary Input (Krstic 2009):
% % PDE         x_{t}     = x_{s1s1} + lam * x;     s1 in [0,1]
% % With BCs    x(t,s1=0) = 0;   
% %             x(t,s1=1) = u(t-tau);
% %
% % Introduce transport equation:
% %     u(t) = x2(t,1+tau),   x2_{t} = x2_{s2}        s2 in [1,1+tau]
% % Then
% %
% % PDE:    x1_{t}         = x1_{s1s1} + lam * x1;
% %         x2_{t}         = x2_{s2};
% % BCs:    x1(t,s1=0)     = 0;
% %         x1(t,s1=1)     = x2(t,s2=1);
% %         x2(t,s2=1+tau) = u(t);
% %
% % Parameters tau and lam can be set.
% % With u = 0 the system is exponentially stable for lam < pi^2.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input
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
% Initial coding DJ - 08/20/2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: second output st, the
%   recommended settings of the stability test; stability_form = 'Q' set in
%   the base workspace; help revised (OUTPUT lists no PDE_b, which is not
%   returned; stability statement added), license added.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2 t

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/21/2026
%evalin('base','stability_dual = 1;');

% Default parameters, overridden by the assignment strings in params
ne = 1; tau = 1;  lam = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
% false: keep the delayed input u(t-tau); true: write the delay as the
% transport state x2
expand_delays = false;
if ~expand_delays
    %%% pde_var input format with delay
    clear stateNameGenerator
    pde_var x control u
    x.vars = s1;        x.dom = [0,1];
    PDE_t = [diff(x,t)==diff(x,s1,2)+lam*x;
             subs(x,s1,0)==0;
             subs(x,s1,1)==subs(u,t,t-tau)];


    % %%% Term-based input format with delay.
    % PDE_t.x{1}.vars = s1;       PDE_t.x{1}.dom = [0,1];
    % PDE_t.u{1}.vars = [];   % Input does not vary in space.
    % 
    % % PDE: x1_{t} = x1_{s1s1} + lam * x1;
    % PDE_t.x{1}.term{1}.x = [1;1];
    % PDE_t.x{1}.term{1}.D = [2; 0];
    % PDE_t.x{1}.term{1}.C = [eye(ne), lam*eye(ne)];
    % 
    % % BC1: 0 = x1(t,s1=0)
    % PDE_t.BC{1}.term{1}.x = 1;
    % PDE_t.BC{1}.term{1}.loc = 0;
    % % BC2: 0 = x1(t,s1=1) - u(t-tau)
    % PDE_t.BC{2}.term{1}.x = 1;          PDE_t.BC{2}.term{2}.u = 1;
    % PDE_t.BC{2}.term{1}.loc = 1;        PDE_t.BC{2}.term{2}.delay = tau;
    %                                     PDE_t.BC{2}.term{2}.C = -1;
    
else
    %%% pde_var input format
    clear stateNameGenerator
    pde_var x1 x2 control u
    x1.vars = s1;       x1.dom = [0,1];
    x2.vars = s2;       x2.dom = [1,1+tau];
    PDE_t = [diff(x1,t)==diff(x1,s1,2)+lam*x1;
             diff(x2,t)==diff(x2,s2);
             subs(x1,s1,0)==0;  subs(x1,s1,1)==subs(x2,s2,1);
             subs(x2,s2,1+tau)==u];


    % %%% Term-based input format
    % % Initialize the state variables and input.
    % PDE_t.x{1}.vars = s1;       PDE_t.x{1}.dom = [0,1];
    % PDE_t.x{2}.vars = s2;       PDE_t.x{2}.dom = [1,1+tau];
    % PDE_t.u{1}.vars = [];   % Input does not vary in space.
    % 
    % % PDE: x1_{t} = x1_{s1s1} + lam * x1;
    % PDE_t.x{1}.term{1}.x = [1;1];
    % PDE_t.x{1}.term{1}.D = [2; 0];
    % PDE_t.x{1}.term{1}.C = [eye(ne), lam*eye(ne)];
    % 
    % % PDE: x2_{t} = x2_{s2}
    % PDE_t.x{2}.term{1}.x = 2;
    % PDE_t.x{2}.term{1}.D = 1;
    % 
    % % BC1: 0 = x1(s1=0)
    % PDE_t.BC{1}.term{1}.x = 1;
    % PDE_t.BC{1}.term{1}.loc = 0;
    % % BC2: 0 = x1(s1=1) - x2(s2=1)
    % PDE_t.BC{2}.term{1}.x = 1;          PDE_t.BC{2}.term{2}.x = 2;
    % PDE_t.BC{2}.term{1}.loc = 1;        PDE_t.BC{2}.term{2}.loc = 1;
    %                                     PDE_t.BC{2}.term{2}.C = -1;
    % % BC3: 0 = x2(s2=1+D) - u;
    % PDE_t.BC{3}.term{1}.x = 2;          PDE_t.BC{3}.term{2}.u = 1;
    % PDE_t.BC{3}.term{1}.loc = 1+tau;      PDE_t.BC{3}.term{2}.C = -1;
end

% Example-specific settings: the light tier                                 % MMP, 11/21/2026
st = lpisettings('light');                                                  % MMP, 11/21/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{krstic2009control,
%   title={Control of an unstable reaction--diffusion PDE with long input delay},
%   author={Krstic, Miroslav},
%   journal={Systems \& Control Letters},
%   volume={58},
%   number={10-11},
%   pages={773--782},
%   year={2009},
%   publisher={Elsevier}
% }