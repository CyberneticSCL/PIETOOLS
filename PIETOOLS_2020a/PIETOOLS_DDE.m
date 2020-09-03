%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_DDE.m     PIETOOLS 2020a
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
close all; clear all; path(pathdef); clc; pvar s theta;
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
% Inputs: A0   - a nx x nx real matrix (required) 
%         tau(i) - These can be an arbitrary sequence of positive increasing
%         numbers. (required)
%
%         A{i} - matrices of dimension n x n (n is arbitrary)
%
%         Ad{i} - matrix (n x n) of polynomial entries in independant variable s (a pvar) 
%
%         B1 - a matrix of dimension n x m (m is arbitrary)
%         B2 - a matrix of dimension n x p (p is arbitrary)
%
%         C1, C1{i} - matrices of dimension q x n (q is arbitrary)
%         C2, C2{i} - matrices of dimension r x n (r is arbitrary)
%
%         D11 - a matrix of dimension q x m
%         D12 - a matrix of dimension q x p
%         D21 - a matrix of dimension r x m
%         D22 - a matrix of dimension r x p
%
%         tau(i) - These can be an arbitrary sequence of positive increasing
%         numbers.
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

% NOTE!!!: By default, PIETOOLS will convert the DDE to a minimum dimension
% DDF representation. This often result in greatly reduced computational
% complexity and for most purposes, this should be left on. However, if you
% wish to use the standard PIE representation of a DDE, then uncomment the
% following line:
%
 %DDE_minimal_rep=0; 

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
% END: USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
eppos2 = 0*1e-4;   % Positivity of Lyapunov Function with respect to spatially distributed states
epneg = 0*1e-5;    % Negativity of Derivative of Lyapunov Function 

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

%settings_PIETOOLS_veryheavy
%settings_PIETOOLS_heavy
settings_PIETOOLS_light
%settings_PIETOOLS_stripped
%settings_PIETOOLS_extreme
 sosineq_on=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END: Linear PI Inequality (LPI) solver settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Convert your DDE system to a PIE
initialize_PIETOOLS_DDE % error checking an preprocessing

if DDE_minimal_rep==1
   % % Conversion from DDE to DDF using a minimal DDF representation.
    minimize_PIETOOLS_DDE
else
   % % Conversion from DDE to DDF using the formulation in our Automatica Paper.
    convert_PIETOOLS_DDE;
end

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


