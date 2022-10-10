%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_DDF.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Started: This file is meant to be run from the Matlab editor 
% window and combines many of the scripts and functions in 
% PIETOOLS for the purpose of studying Differential-Difference Equations (DDFs).
%
%
% USER INPUT: The first set of commands inputs a DDF using the format definition in the header.
% This user input can be replaced with a call to the library file
% examples_DDF_library_PIETOOLS which are all in the correct format. At the
% moment, however, the list of examples in this file is relatively sparse.
% For most people, you should use PIETOOLS_DDE.
%
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
% addpath(genpath('.')) % makes sure any local version  of SOSTOOLS in the
% current folder is at the head of the path (Troubleshooting purposes only)
% close all; clear all; path(pathdef); clc; 
clear;
pvar s theta; sosineq_on=0;
stability=0; stability_dual=0; Hinf_gain=0; Hinf_gain_dual=0; Hinf_control=0; Hinf_estimator=0; 
DDE2DDF=0;DDF_min=1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN: USER INPUT: Define the DDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main executable file for studying DDF systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program determines stability, hinf-norm and designs hinf-optimal
% observer of a linear DDF which is defined in the format given below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5-Equation Universal Formulation of the DDF:
%
%	\dot{x}(t)=A0x(t)+B1w(t)+B2u(t)+Bv v(t)
% where w(t) is the vector of disturbances, u(t) is the vector of control
% inputs, and v(t) is the output from the delayed channels
%
% r_i are the output signals to be delayed - dimension rdim{i}
% r_i(t)=Cr{i}x(t)+Br1{i}w(t)+Br2{i}u(t)+Drv{i}v(t)
%
% v is the output of the delayed channels
% v(t)=\sum_{i=1}^K Cv{i} r_i(t-\tau_i)+\sum_{i=1}^K \int_{-\tau_i}^0Cvd{i}(s) r_i(t+s)ds
%
%
% The regulated output is of the form
%	z(t)=C1 x(t)+D11 w(t)+D12 u(t)+D1v v(t)
%
% The sensed output is of the form
%	y(t)=C2 x(t)+D21 w(t)+D22 u(t)+D2v v(t)
%
% u(t) is the controller input and is of the form
% u(t)=K0 x(t) + \sum_i K1{i}x(t-tau(i)) + \sum_i \int_{-\tau(i)}^0 K2{i}(s)x(t-s)ds
% where A0, Ai{i}, Bi, Cij{i}, Dij and tau(i) are all user inputs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-Defined Inputs for specifying the dynamics:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Inputs: 
% These are grouped by the 5 equations in the DDF formulation
% which they are used to represent. All elements should be stored to a DDF data structure. The
% possible elements of this data structure which can be defined are as
% follows.
%
%         DDF.tau(i) - These can be an arbitrary sequence of positive increasing
%         numbers. (required)
%
%         DDF.A0 - a matrix of dimension n x n (optional)
%         DDF.Ai{i} - a cell of ntau(# of delays) matrices of dimension n x n (optional)
%         DDF.B1 - a matrix of dimension n x nw (optional)
%         DDF.B2 - a matrix of dimension n x nu (optional)
%         DDF.Bv - a matrix of dimension n x nv (optional)
%
%         DDF.Cr{i} - a cell of ntau(# of delays) matrices of dimension nr{i} x n (optional)
%         DDF.Br1{i} - a cell of ntau(# of delays) matrices of dimension nr{i} x nw (optional)
%         DDF.Br2{i} - a cell of ntau(# of delays) matrices of dimension nr{i} x nu (optional)
%         DDF.Drv{i} - a cell of ntau(# of delays) matrices of dimension nr{i} x nv (optional)
%
%         DDF.Cv{i} - a cell of ntau(# of delays) matrices of dimension nv x nr{i} (optional)
%         DDF.Cdv{i} - a cell of ntau(# of delays) matrix-valued polynomials of dimension nv x nr{i} (optional)
%
%         DDF.C1 - a matrix of dimension nz x n (optional)
%         DDF.D11 - a matrix of dimension nz x nw (optional)
%         DDF.D12 - a matrix of dimension nz x nu (optional)
%         DDF.D1v - a matrix of dimension nz x nv (optional)
%
%         DDF.C2 - a matrix of dimension ny x n (optional)
%         DDF.D21 - a matrix of dimension ny x nw (optional)
%         DDF.D22 - a matrix of dimension ny x nu (optional)
%         DDF.D2v - a matrix of dimension ny x nv (optional)
%
%         The selection of tools which are desired: set the desired options to 1.
% stability=1
% stability_dual=1
% Hinf_gain=1
% Hinf_gain_dual=1
% Hinf_control=1
% Hinf_estimator=1
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Library Call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To test PIETOOLS on one of the standard example problems, consult the
% library file examples_DDF_library_PIETOOLS.m and uncomment the example
% you wish to test. Otherwise, comment the following line and enter your
% own DDF.

examples_NDSDDF_library_PIETOOLS



%DDF_min=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using DDE examples as DDFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The library of DDF examples is still relatively sparse.
% As an alternative to entering a custom DDF, the user may enter any of the
% DDEs in the library examples_DDE_library_PIETOOLS and use the option
% DDE2DDF=1 to convert the DDE to an equivalent DDF, which is then
% converted to a PIE.
% NOTE 1: Converting a DDE into a DDF using DDE2DDF=1 will lose any special
% structure of the system which might otherwise have been exploited to
% reduce computational complexity. 
% NOTE 2: As an alternative, we have a optimized DDE to DDF converter which
% identifies any special structure in the DDE and exploits this to find the
% minimum dimension DDF representation. In fact, the default in the DDE
% executable is to identify this structure and convert to a DDF and from
% there to a PIE. If you have specified a DDE and wish to exploit this
% structure, use the following option.

% examples_DDE_library_PIETOOLS
% DDE2DDF=1;

%%%%%%%%%%%%% Choose Problem to Solve %%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have not already chosen the problem you want to solve, then
% uncomment one of the following lines:
%

% stability=0;             % Test Lyapunov Stability
% stability_dual=0;        % Test Lyapunov Stability using an alternative Duality theorem
% Hinf_gain=0;             % Find a minimal bound on the Hing gain of the system
% Hinf_gain_dual=0;        % Find a minimal bound on the Hing gain of the system using an alternative duality theorem
% Hinf_control=0;          % Find a minimal bound on the Hinf optimal control problem
% Hinf_estimator=0;        % Find a minimal bound on the Hinf optimal observer problem


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

sosineq_on=0; % binary variable indicating whether to use ineqaulity or equality constraint

% settings_PIETOOLS_veryheavy
% settings_PIETOOLS_heavy
 settings_PIETOOLS_light;
% settings_PIETOOLS_stripped
% settings_PIETOOLS_extreme
% settings_PIETOOLS_custom


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
settings.epneg = 0*1e-5;    % Negativity of Derivative of Lyapunov Function %%% Default Setttings: %%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END: Linear PI Inequality (LPI) solver settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Convert your DDF system to a PIE
if exist('NDS','var')
     NDS=initialize_PIETOOLS_NDS(NDS);
     DDF=convert_PIETOOLS_NDS2DDF(NDS);
end
if DDE2DDF==1
    if DDF_min==0
        DDE=initialize_PIETOOLS_DDE(DDE);
        PIE=convert_PIETOOLS_DDE2PIE(DDE);
    else
        DDE=initialize_PIETOOLS_DDE(DDE);
        DDE=minimize_PIETOOLS_DDE2DDF(DDE);
        PIE=convert_PIETOOLS_DDE2PIE(DDE);
    end
else
    if DDF_min==1
        % % Construct a minimal realizatio of the DDF and then convert to a PIE
        DDF = initialize_PIETOOLS_DDF(DDF);
        DDF = minimize_PIETOOLS_DDF(DDF);
        PIE = convert_PIETOOLS_DDF(DDF,'pie');
    else
        DDF = initialize_PIETOOLS_DDF(DDF);
        PIE = convert_PIETOOLS_DDF(DDF,'pie');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Calling the Executives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    settings.options1.sep = 1;  % Pop needs to separable in order to construct controller - SS,MP
    settings.options12.sep = 1; 

    [prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);
end
if Hinf_estimator==1
    settings.options1.sep = 1;  % Pop needs to separable in order to construct controller - SS,MP
    settings.options12.sep = 1; 
    [prog, L, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT: Format for defining an Neutral Delay System (NDS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This root script may also be used to define and analyse an NDS given in the format below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% \dot{x}(t)=  A0x(t) + \sum_i Ai{i}  x(t-tau(i))     + \sum_i \int_{-tau(i}^0 Adi{i}(s)  x(t+s)ds
%                       \sum_i Ei{i}  \dot x(t-tau(i))+ \sum_i \int_{-tau(i}^0 Edi{i}(s)  \dot x(t+s)ds
%            + B1w(t) + \sum_i B1i{i} w(t-tau(i))     + \sum_i \int_{-tau(i}^0 B1di{i}(s) w(t+s)ds
%            + B2u(t) + \sum_i B2i{i} u(t-tau(i))     + \sum_i \int_{-tau(i}^0 B2di{i}(s) u(t+s)ds
% where w(t) is the vector of disturbances
%
% The regulated output is of the form
% z(t) =  C1x(t)  + \sum_i C1i{i}  x(t-tau(i))      + \sum_i \int_{-tau(i}}^0 C1d{i}(s)   x(t+s)ds
%                 + \sum_i E1i{i}  \dot x(t-tau(i)) + \sum_i \int_{-tau(i}^0 E1di{i}(s)   \dot x(t+s)ds
%       + D11w(t) + \sum_i D11i{i} w(t-tau(i))      + \sum_i \int_{-tau(i)}^0 D11di{i}(s) w(t+s)ds
%       + D12u(t) + \sum_i D12i{i} u(t-tau(i))      + \sum_i \int_{-tau(i)}^0 D12di{i}(s) u(t+s)ds
%
% The sensed output is of the form
% y(t) =  C2x(t)  + \sum_i C2i{i}  x(t-tau(i))      + \sum_i \int_{-tau(i}}^0 C2d{i}(s)   x(t+s)ds
%                 + \sum_i E2i{i}  \dot x(t-tau(i)) + \sum_i \int_{-tau(i}^0 E2di{i}(s)   \dot x(t+s)ds
%       + D21w(t) + \sum_i D21i{i} w(t-tau(i))      + \sum_i \int_{-tau(i)}^0 D21di{i}(s) w(t+s)ds
%       + D22u(t) + \sum_i D22i{i} u(t-tau(i))      + \sum_i \int_{-tau(i)}^0 D22di{i}(s) u(t+s)ds
%
% where A0, Ai{i}, Adi{i}, Ei{i}, E1i{i}, E2i{i}, Edi{i}, E1di{i}, E2di{i}, Bi{i}, B1i{i}, B2i{i}, B1di{i}, B2di{i}, tau(i), etc. are all user inputs.
%
% Inputs: All elements should be stored to a NDS data structure. The
% possible elements of this data structure which can be defined are as
% follows.
%
%         NDS.A0   - a nx x nx real matrix (required) 
%         NDS.tau(i) - These can be an arbitrary sequence of positive increasing
%         numbers. (required)
%
%         NDS.A{i} - matrices of dimension n x n (n is arbitrary)
%         NDS.Ad{i} - matrix (n x n) of polynomial entries in independant variable s (a pvar) 
%
%         NDS.Ei{i}  - matrices of dimension n x n (n is arbitrary)
%         NDS.Edi{i} - matrix function of variable pvar s of dimension n x n (n is arbitrary)
%
%         NDS.B1, NDS.B1i{i} - a matrix of dimension n x m (m is arbitrary)
%         NDS.B2, NDS.B2i{i} - a matrix of dimension n x p (p is arbitrary)
%         NDS.B1di{i} - matrix function of variable pvar s of dimension n x m (q is arbitrary)
%         NDS.B2di{i} - matrix function of variable pvar s of dimension n x p (r is arbitrary)
%
%         NDS.C1, NDS.C1i{i} - matrices of dimension q x n (q is arbitrary)
%         NDS.C2, NDS.C2i{i} - matrices of dimension r x n (r is arbitrary)
%         NDS.C1di{i} - matrix function of variable pvar s of dimension q x n (q is arbitrary)
%         NDS.C2di{i} - matrix function of variable pvar s of dimension r x n (r is arbitrary)
%
%         NDS.D11, NDS.D11i{i} - a matrix of dimension q x m
%         NDS.D12, NDS.D12i{i} - a matrix of dimension q x p
%         NDS.D21, NDS.D21i{i} - a matrix of dimension r x m
%         NDS.D22, NDS.D22i{i} - a matrix of dimension r x p
%         NDS.D11di{i} - matrix function of variable pvar s of dimension q x m
%         NDS.D12di{i} - matrix function of variable pvar s of dimension q x p
%         NDS.D21di{i} - matrix function of variable pvar s of dimension r x m
%         NDS.D22di{i} - matrix function of variable pvar s of dimension r x p

