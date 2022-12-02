%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_DDE.m     PIETOOLS 2021b
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
%  LPI SOLVER SETTINGS: The first set of commands in the header specify default settings for the
%  Matlab environment and the LPI solver. This includes a call to one of 5
%  settings scripts, and a declaration for how strictly the Lyapunov
%  inequalities are enforced. If you are not solving an LPI, you may ignore
%  all of these commands.
%
% The third set of commands converts the DDE to a PIE either directly or by
% converting of a minimal Delay-Differential Equation and hence to a PIE.
%
% The fourth and final set of commands uses a script which sets up and solves an LPI as
% determined by the LPI settings and 
%
close all; clear all; %path(pathdef); 
clc; pvar s theta;
stability=0; stability_dual=0; Hinf_gain=0; Hinf_gain_dual=0; Hinf_control=0; Hinf_estimator=0;
sosineq_on=1; DDE_minimal_rep=1;
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
%         DDE.A{i} - matrices of dimension n x n (n is arbitrary)
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
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Library Call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To test PIETOOLS on one of the standard example problems, consult the
% library file examples_dde_library_PIETOOLS.m and uncomment the example
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
 uinput.ic.ODE=[5,5,5];
 solution=PIESIM(DDE,opts,uinput);
 end
 
%------------------------------------------------
% Description of solution field - output of the simulations
%------------------------------------------------
% solution.tf - final time of the solution 
% solution.final.pde - distributed state solution at a final time tf: matrix of
% dimension (N+1) x ns, ns - total number of distributed states
% solution.final.ode - ode solution at a final time tf: array of dimensions
% nx x 1, nx - total numer of ODE states
% solution.timedep.ode - array of size nx x Nsteps - time-dependent solution of nx ODE states
% solution.timedep.pde - array of size (N+1) x ns x Nsteps - time-dependent
% solution of ns distributed states
% solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
%---------------------------------------------------

    
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
% Uncomment your desired settings script:
%
% settings_PIETOOLS_veryheavy is for difficult problems which require complex
% Lyapunov/storage functions.
%
% settings_PIETOOLS_heavy is reliable/accurate and should work for almost 
% all problems - at some cost of computational complexity
%
% settings_PIETOOLS_light works well for most problems, but but may
% sacrifice accuracy in the second decimal point
%
% settings_PIETOOLS_stripped sacrifices accuracy
%
% settings_PIETOOLS_extreme is often infeasible


settings = lpisettings('light');
settings.sosineq_on=0; % binary variable indicating whether to use ineqaulity or equality constraint


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Default Setttings: %%%%%%
settings.sos_opts.solver='sedumi';
% Other optional SDP solvers supported by PIETOOLS
% settings.sos_opts.solver ='mosek';
% settings.sos_opts.solver='sdpnalplus';
% settings.sos_opts.solver='sdpt3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strictness of Inequalities:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

settings.eppos = 1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
settings.eppos2 = 0*1e-4;   % Positivity of Lyapunov Function with respect to spatially distributed states
settings.epneg = 0*1e-5;    % Negativity of Derivative of Lyapunov Function 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END: Linear PI Inequality (LPI) solver settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Convert your DDE system to a PIE
DDE=initialize_PIETOOLS_DDE(DDE); % error checking an preprocessing

if DDE_minimal_rep==1
   % % Conversion from DDE to DDF using a minimal DDF representation.
    DDE = initialize_PIETOOLS_DDE(DDE); % error checking and preprocessing
    DDF = minimize_PIETOOLS_DDE2DDF(DDE);
    PIE = convert_PIETOOLS_DDF(DDF,'pie');  
else
   % % Conversion from DDE to DDF using the formulation in our Automatica Paper.
    PIE=convert_PIETOOLS_DDE2PIE(DDE);
end


%%% Construct and solve the LPI corresponding to the problem you specified earlier
if stability==1
    [prog, P] = PIETOOLS_stability(PIE,settings);
end
if stability_dual==1
    [prog, P] = PIETOOLS_stability_dual(PIE,settings);
end
if Hinf_gain==1
    [prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);
end
if Hinf_gain_dual==1
    [prog, P, gamma] = PIETOOLS_Hinf_gain_dual(PIE,settings);
end
if Hinf_control==1
    [prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);
end
if Hinf_estimator==1
    [prog, L, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings);
end


