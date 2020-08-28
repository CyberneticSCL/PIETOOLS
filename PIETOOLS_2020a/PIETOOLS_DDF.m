%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_DDF.m     PIETOOLS 2020a
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
close all; clear all; path(pathdef); clc; pvar s theta; sosineq_on=0;
stability=0; stability_dual=0; Hinf_gain=0; Hinf_gain_dual=0; Hinf_control=0; Hinf_estimator=0; 
DDE2DDF_min=0;
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
% User Inputs: These are grouped by the 5 equations in the DDF formulation
% which they are used to represent
%
%         tau(i) - These can be an arbitrary sequence of positive increasing
%         numbers. (required)
%
%         A0 - a matrix of dimension n x n (optional)
%         Ai{i} - a cell of ntau(# of delays) matrices of dimension n x n (optional)
%         B1 - a matrix of dimension n x nw (optional)
%         B2 - a matrix of dimension n x nu (optional)
%         Bv - a matrix of dimension n x nv (optional)
%
%         Cr{i} - a cell of ntau(# of delays) matrices of dimension nr{i} x n (optional)
%         Br1{i} - a cell of ntau(# of delays) matrices of dimension nr{i} x nw (optional)
%         Br2{i} - a cell of ntau(# of delays) matrices of dimension nr{i} x nu (optional)
%         Drv{i} - a cell of ntau(# of delays) matrices of dimension nr{i} x nv (optional)
%
%         Cv{i} - a cell of ntau(# of delays) matrices of dimension nv x nr{i} (optional)
%         Cdv{i} - a cell of ntau(# of delays) matrix-valued polynomials of dimension nv x nr{i} (optional)
%
%         C1 - a matrix of dimension nz x n (optional)
%         D11 - a matrix of dimension nz x nw (optional)
%         D12 - a matrix of dimension nz x nu (optional)
%         D1v - a matrix of dimension nz x nv (optional)
%
%         C2 - a matrix of dimension ny x n (optional)
%         D21 - a matrix of dimension ny x nw (optional)
%         D22 - a matrix of dimension ny x nu (optional)
%         D2v - a matrix of dimension ny x nv (optional)
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

examples_DDF_library_PIETOOLS

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

%examples_DDE_library_PIETOOLS
%DDE2DDF_min=1;

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
%%% Default Setttings: %%%%%%
sos_opts.solver='sedumi';
% Other optional SDP solvers supported by PIETOOLS
% sos_opts.solver ='mosek';
% sos_opts.solver='sdpnalplus';
% sos_opts.solver='sdpt3';

% sosineq_on=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strictness of Inequalities:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eppos = 1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
eppos2 = 1*1e-6;   % Positivity of Lyapunov Function with respect to spatially distributed states
epneg = 0*1e-5;    % Negativity of Derivative of Lyapunov Function in both current and distributed state -  >0 if exponential stability desired

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational Complexity Toggles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: DDEs and DDFs do not require high degrees and hence the
% settings_PIETOOLS_light script is adequate for most applications

% settings_PIETOOLS_veryheavy
% settings_PIETOOLS_heavy
  settings_PIETOOLS_light
% settings_PIETOOLS_stripped
% settings_PIETOOLS_extreme
% settings_PIETOOLS_custom


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END: Linear PI Inequality (LPI) solver settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Convert your ODE-PDE system to a PIE

if DDE2DDF_min==1
% % Conversion from DDE to DDF using a minimal DDF representation
    initialize_PIETOOLS_DDE % error checking and preprocessing
    minimize_PIETOOLS_DDE
else
    initialize_PIETOOLS_DDF;
    convert_PIETOOLS_DDF;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Calling the Executives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if stability==1
    executive_PIETOOLS_stability;
end
if stability_dual==1
    executive_PIETOOLS_stability_dual;
end
if Hinf_gain==1
    executive_PIETOOLS_Hinf_gain;
end
if Hinf_gain_dual==1
    executive_PIETOOLS_Hinf_gain_dual;
end
if Hinf_control==1
    executive_PIETOOLS_Hinf_control;
end
if Hinf_estimator==1
    executive_PIETOOLS_Hinf_estimator;
end

