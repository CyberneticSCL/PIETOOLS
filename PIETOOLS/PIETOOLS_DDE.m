%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_DDE.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Getting Started: This file is meant to be run from the Matlab editor 
% window and combines many of the scripts and functions in 
% PIETOOLS for the purpose of studying Delay-Differential Equations (DDEs).
%
%
% USER INPUT: The first set of commands inputs a DDE using the format definition in the header.
% This user input can be replaced with a call to the library file
% examples_DDE_library_PIETOOLS which are all in the correct format.
%
% With simulate = 'on', the DDE is then simulated by PIESIM.
%
%  LPI SOLVER SETTINGS: The first set of commands in the header specify default settings for the
%  Matlab environment and the LPI solver. This includes the settings of a
%  tier, lpisettings('light'), replaced by the settings st of the
%  library example when it sets them, and a declaration for how strictly the Lyapunov
%  inequalities are enforced. If you are not solving an LPI, you may ignore
%  all of these commands.
%
% The third set of commands converts the DDE to a PIE either directly or
% through its minimal Differential-Difference Equation (DDF) representation
% (DDE_minimal_rep = 1, the default).
%
% The fourth and final set of commands runs the executives that the flags
% select (stability, stability_dual, Hinf_gain, Hinf_gain_dual,
% Hinf_control, Hinf_estimator, H2_norm, H2_norm_dual, H2_control,
% H2_estimator), each of which sets up and solves an LPI as determined by
% the LPI settings and, for stability and the H-infinity gain, by the forms
% stability_form and Hinf_gain_form.
%
% INPUT
% - none; the DDE (structure DDE, in the format below), the flags of the
%   problems to solve and, for some examples, the LPI settings st and the
%   forms stability_form and Hinf_gain_form come from the library call, or
%   from declarations in its place.
%
% OUTPUT
% - none; the workspace receives the DDE, the simulation output solution,
%   the PIE, the LPI settings, and the outputs of the executives run (prog,
%   P and, per problem, gamma, K, L, Z, W, R, Q).
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_DDE
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
% Initial coding PIETOOLS Team (PIETOOLS 2021b; in the repository since 12/31/2021)
% YP, DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPIs are solved
%   by the container executives, with the settings of the 'light' tier or
%   the settings st of a library example that sets them; stability_form =
%   'P' and Hinf_gain_form = 'Q' set the forms of the stability and
%   H-infinity gain tests (the dual appends 'd'); the DDE is initialized
%   once; the direct conversion calls convert_PIETOOLS_DDE2PIE_legacy; the
%   PIESIM initial condition is uinput.ic; help and comments revised.
% MMP, 11/27/2026: copvar, cdopvar renamed opvar, dopvar, with the
%   functions of the rename table; names only.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
%
close all; clear; %path(pathdef); 
clc;
stability=0; stability_dual=0; 
Hinf_gain=0; Hinf_gain_dual=0; H2_norm=0; H2_norm_dual=0;
Hinf_control=0; Hinf_estimator=0; H2_control=0; H2_estimator=0;
stability_form='P';     % storage form; a library block may set 'Q'         % YP, DB, MMP, 11/21/2026
Hinf_gain_form='Q';     % gain form; a library block may set 'P'            % YP, DB, MMP, 11/21/2026

DDE_minimal_rep=1;
% sosineq_on=1; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN: USER INPUT: Define the DDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main executable file for studying DDE systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program determines stability, hinf-norm and designs hinf-optimal
% observers and controllers of a DDE which is defined in the format given below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% \dot{x}(t)=  A0x(t) + \sum_i Ai{i}  x(t-tau(i)) + \sum_i \int_{-tau(i}^0 Adi{i}(s)x(t+s)ds
%            + B1w(t) + \sum_i B1i{i} w(t-tau(i)) + \sum_i \int_{-tau(i}^0 B1di{i}(s) w(t+s)ds
%            + B2u(t) + \sum_i B2i{i} u(t-tau(i)) + \sum_i \int_{-tau(i}^0 B2di{i}(s) u(t+s)ds
% where w(t) is the vector of disturbances
%
% The regulated output is of the form
% z(t) =  C1x(t) + \sum_i C1i{i}x(t-tau(i))   + \sum_i \int_{-tau(i}}^0 C1d{i}(s) x(t+s)ds
%       + D11w(t) + \sum_i D11i{i} w(t-tau(i)) + \sum_i \int_{-tau(i)}^0 D11di{i}(s) w(t+s)ds
%       + D12u(t) + \sum_i D12i{i} u(t-tau(i)) + \sum_i \int_{-tau(i)}^0 D12di{i}(s) u(t+s)ds
%
% The sensed output is of the form
% y(t) =  C2x(t) + \sum_i C2i{i}x(t-tau(i))   + \sum_i \int_{-tau(i}}^0 C2d{i}(s) x(t+s)ds
%       + D21w(t) + \sum_i D21i{i} w(t-tau(i)) + \sum_i \int_{-tau(i)}^0 D21di{i}(s) w(t+s)ds
%       + D22u(t) + \sum_i D22i{i} u(t-tau(i)) + \sum_i \int_{-tau(i)}^0 D22di{i}(s) u(t+s)ds
%
% where A0, Ai{i}, Bi, Cij{i}, Dij and tau(i) are all user inputs.
%
% Inputs: All elements should be stored to a DDE data structure. The
% possible elements of this data structure which can be defined are as
% follows.
%%%%%%%%%%%%%%%%%% Required Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         DDE.A0   - a nx x nx real matrix (required) 
%         DDE.tau(i) - These can be an arbitrary sequence of positive increasing
%         numbers. (required)
%
%%%%%%%%%%%%%%%%%%%%%%%% Optional Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         DDE.Ai{i} - matrices of dimension n x n (n is arbitrary)
%         DDE.Ad{i} - matrix (n x n) of polynomial entries in independant variable s (a pvar) 
%
%         DDE.B1, DDE.B1i{i} - a matrix of dimension n x m (m is arbitrary)
%         DDE.B2, DDE.B2i{i} - a matrix of dimension n x p (p is arbitrary)
%         DDE.B1di{i} - matrix function of variable pvar s of dimension n x m (q is arbitrary)
%         DDE.B2di{i} - matrix function of variable pvar s of dimension n x p (r is arbitrary)
%
%         DDE.C1, DDE.C1i{i} - matrices of dimension q x n (q is arbitrary)
%         DDE.C2, DDE.C2i{i} - matrices of dimension r x n (r is arbitrary)
%         DDE.C1di{i} - matrix function of variable pvar s of dimension q x n (q is arbitrary)
%         DDE.C2di{i} - matrix function of variable pvar s of dimension r x n (r is arbitrary)
%
%         DDE.D11, DDE.D11i{i} - a matrix of dimension q x m
%         DDE.D12, DDE.D12i{i} - a matrix of dimension q x p
%         DDE.D21, DDE.D21i{i} - a matrix of dimension r x m
%         DDE.D22, DDE.D22i{i} - a matrix of dimension r x p
%         DDE.D11di{i} - matrix function of variable pvar s of dimension q x m
%         DDE.D12di{i} - matrix function of variable pvar s of dimension q x p
%         DDE.D21di{i} - matrix function of variable pvar s of dimension r x m
%         DDE.D22di{i} - matrix function of variable pvar s of dimension r x p
%
%
%
%
%
%         The selection of tools which are desired: set the desired options to 1.
% stability=1
% stability_dual=1
% Hinf_gain=1
% Hinf_gain_dual=1
% Hinf_control=1
% Hinf_estimator=1
%
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Library Call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To test PIETOOLS on one of the standard example problems, consult the
% library file examples_DDE_library_PIETOOLS.m and uncomment the example
% you wish to test. Otherwise, comment the following line and enter your
% own DDE.

examples_DDE_library_PIETOOLS

simulate='on';

if strcmp(simulate,'on')
%------------------------------------------------
% Set-up default options for simulations
% The following default options are used in simulation
% N=8 - number of discretization points
% dt=1e-3 - time step increment
% tf=10 - final time of the solution
% Options for disturbance shape (if applicable): 'constant', 'sin' or 'sinc' 
% Options for plotting: select 'yes' to plot ODE states, 'no' if not
%-------------------------------------------------
 opts.N=8;
 opts.dt=1e-2;
 opts.tf=10;
 opts.dist='constant';
 opts.plot='yes';
 uinput.ic=[5,5,5];
 solution=PIESIM(DDE,opts,uinput);
 end
 
%------------------------------------------------
% solution: 'struct' with fields tf (the final time), final (the states,
% observed and regulated outputs at tf) and timedep (their values at the
% time steps dtime)
%------------------------------------------------

    
% NOTE!!!: By default, PIETOOLS will convert the DDE to a minimum dimension
% DDF representation. This often result in greatly reduced computational
% complexity and for most purposes, this should be left on. However, if you
% wish to use the standard PIE representation of a DDE, then uncomment the
% following line:
%
% DDE_minimal_rep=0; 

%%%%%%%%%%%%% Choose Problem to Solve %%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have not already chosen the problem you want to solve, then
% uncomment one of the following lines:
%

% stability=1;             % Test Lyapunov Stability
% stability_dual=1;        % Test Lyapunov Stability using an alternative Duality theorem
% Hinf_gain=1;             % Find a minimal bound on the Hing gain of the system
% Hinf_gain_dual=1;        % Find a minimal bound on the Hing gain of the system using an alternative duality theorem
% Hinf_control=1;          % Find a minimal bound on the Hinf optimal control problem
% Hinf_estimator=1;        % Find a minimal bound on the Hinf optimal observer problem


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END: USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN: Linear PI Inequality (LPI) solver settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational Complexity Toggles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose the tier of lpisettings below:
%
% 'veryheavy' is for difficult problems which require complex
% Lyapunov/storage functions.
%
% 'heavy' is reliable/accurate and should work for almost
% all problems - at some cost of computational complexity
%
% 'light' works well for most problems, but but may
% sacrifice accuracy in the second decimal point
%
% 'stripped' sacrifices accuracy
%
% 'extreme' is often infeasible


sosineq_on=0; % binary variable indicating whether to use ineqaulity or equality constraint
% settings = lpisettings('light');                                          % YP, DB, MMP, 11/21/2026 (was)
settings = lpisettings('light');      % tiers of the container executives   % YP, DB, MMP, 11/21/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Default Setttings: %%%%%%
% % OPTIONAL: Uncomment to declare SDP solver to use (defaults to solver on path)
% settings.sos_opts.solver = 'sedumi';
% settings.sos_opts.solver = 'mosek';
% settings.sos_opts.solver = 'sdpnalplus';
% settings.sos_opts.solver = 'sdpt3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strictness of Inequalities:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

settings.eppos = 1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
settings.eppos2 = 0*1e-4;   % Positivity of Lyapunov Function with respect to spatially distributed states
settings.epneg = 0*1e-5;    % Negativity of Derivative of Lyapunov Function 
% The settings st of the library example, if it sets them, replace these    % YP, DB, MMP, 11/21/2026
if exist('st','var'),   settings = st;  end                                 % YP, DB, MMP, 11/21/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END: Linear PI Inequality (LPI) solver settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Convert your DDE system to a PIE
DDE=initialize_PIETOOLS_DDE(DDE); % error checking an preprocessing

if DDE_minimal_rep==1
   % % Conversion from DDE to DDF using a minimal DDF representation.
%    DDE = initialize_PIETOOLS_DDE(DDE); % error checking and preprocessing
    DDF = minimize_PIETOOLS_DDE2DDF(DDE);
    PIE = convert_PIETOOLS_DDF(DDF,'pie');  
else
   % % Conversion from DDE directly to PIE (every state and input delayed) using the formulation in our Automatica Paper.
    PIE=convert_PIETOOLS_DDE2PIE_legacy(DDE);
end


%%% Construct and solve the LPI corresponding to the problem you specified earlier
% The container executives (PIETOOLS_*) take the opvar PIE that both   % YP, DB, MMP, 11/21/2026
% converters return; outputs as named below, then an info struct            % YP, DB, MMP, 11/21/2026
if stability==1
%   [prog, P] = PIETOOLS_PIE2PDEstability(PIE,settings);                    % YP, DB, MMP, 11/21/2026 (was)
    [prog, P] = PIETOOLS_stability(PIE,settings,stability_form);            % YP, DB, MMP, 11/21/2026
end
if stability_dual==1
%   [prog, P] = PIETOOLS_PIE2PDEstability_dual(PIE,settings);               % YP, DB, MMP, 11/21/2026 (was)
    [prog, P] = PIETOOLS_stability(PIE,settings,[stability_form 'd']);      % YP, DB, MMP, 11/21/2026
end
if Hinf_gain==1
%   [prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);                    % YP, DB, MMP, 11/21/2026 (was)
    [prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings,Hinf_gain_form); % YP, DB, MMP, 11/21/2026
end
if Hinf_gain_dual==1
%   [prog, P, gamma] = PIETOOLS_Hinf_gain_dual(PIE,settings);               % YP, DB, MMP, 11/21/2026 (was)
    [prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings,[Hinf_gain_form 'd']); % YP, DB, MMP, 11/21/2026
end
if Hinf_control==1
%   [prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);           % YP, DB, MMP, 11/21/2026 (was)
    [prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);           % YP, DB, MMP, 11/21/2026
end
if Hinf_estimator==1
%   [prog, L, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings);         % YP, DB, MMP, 11/21/2026 (was)
    [prog, L, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings);         % YP, DB, MMP, 11/21/2026
end
if H2_norm==1
%   [prog, W, gamma, R, Q] = PIETOOLS_H2_norm_o(PIE,settings);              % YP, DB, MMP, 11/21/2026 (was)
    [prog, W, gamma, R, Q] = PIETOOLS_H2_norm(PIE,settings,'o');            % YP, DB, MMP, 11/21/2026
end
if H2_norm_dual==1
%   [prog, W, gamma, R, Q] = PIETOOLS_H2_norm_c(PIE,settings);              % YP, DB, MMP, 11/21/2026 (was)
    [prog, W, gamma, R, Q] = PIETOOLS_H2_norm(PIE,settings,'c');            % YP, DB, MMP, 11/21/2026
end
if H2_control==1
%   [prog, K, gamma, P, Z] = PIETOOLS_H2_controller(PIE,settings);          % YP, DB, MMP, 11/21/2026 (was)
    [prog, K, gamma, P, Z] = PIETOOLS_H2_control(PIE,settings);             % YP, DB, MMP, 11/21/2026
end
if H2_estimator==1
%   [prog, L, gamma, P, Z] = PIETOOLS_H2_estimator(PIE,settings);           % YP, DB, MMP, 11/21/2026 (was)
    [prog, L, gamma, P, Z] = PIETOOLS_H2_estimator(PIE,settings);           % YP, DB, MMP, 11/21/2026
end



