%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_PDE.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Started: This file is meant to be run from the Matlab editor
% window and combines many of the scripts and functions in
% PIETOOLS for the purpose of studying Partial-Differential Equations (PDEs).
%
%
% USER INPUT: The first set of commands inputs a PDE using the format definition in the header.
% This user input can be replaced with a call to the library file
% examples_PDE_library_PIETOOLS which are all in the correct format.
%
%  LPI SOLVER SETTINGS: The first set of commands in the header specify default settings for the
%  Matlab environment and the LPI solver. This includes a call to one of 5
%  settings scripts, and a declaration for how strictly the Lyapunov
%  inequalities are enforced. If you are not solving an LPI, you may ignore
%  all of these commands.
%
% The third set of commands converts the PDE to a PIE.
%
% The fourth and final set of commands uses a script which sets up and solves an LPI as
% determined by the LPI settings and
%
% addpath(genpath('.')) % makes sure any local version  of SOSTOOLS in the
% current folder is at the head of the path (Troubleshooting purposes only)
close all; clear all; path(pathdef); clc; pvar s theta; sosineq_on=0;
stability=0; stability_dual=0; Hinf_gain=0; Hinf_gain_dual=0; Hinf_control=0; Hinf_estimator=0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN: USER INPUT: Define the PDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main executable file for studying coupled ODE-PDE systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program determines stability, hinf-norm and designs hinf-optimal
% observer of a linear coupled ODE-PDE which is defined in the format given below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Definition:
% \dot [xo(t) ]  =A xo(t)+ E0 [x_2(a) + int(Ea(s) [x_1(t,s)] ds,a,b) + int(Eb(s) [x_2s(t,s)] ds,a,b) + int(Ec(s)[x_2ss(t,s)] ds,a,b) + [B11] w(t) + [B12] u(t)
%                              x_2(b)             [x_2(t,s)]                     [x_3s(t,s)]
%                              x_3(a)             [x_3(t,s)]
%                              x_3(b)
%                              x_3s(a)
%                              x_3s(b)]
%
% \dot [x(t,s)] = E(s) xo(t)+ A0(s) [x_1(t,s)] + A1(s) [x_2s(t,s)] + A2(s) [x_3ss(t,s)]+ [B21(s)]w(t)+ [B22(s)]u(t)
%                                   [x_2(t,s)]         [x_3s(t,s)]
%                                   [x_3(t,s)]
%
% z(t) = C1 xo(t) + C10[x_2(a) + int(Ca1(s)[x_1(t,s)] ds,a,b) + int(Cb1(s)[x_2s(t,s)] ds,a,b)+ int(Cc1(s)[x_2ss(t,s)] ds,a,b) + [D11]w(t) + [D12]u(t)
%                       x_2(b)             [x_2(t,s)]                     [x_3s(t,s)]
%                       x_3(a)             [x_3(t,s)]
%                       x_3(b)
%                       x_3s(a)
%                       x_3s(b)]
%
% y(t) = C2 xo(t) + C20[x_2(a) + int(Ca2(s)[x_1(t,s)] ds,a,b) + int(Cb2(s)[x_2s(t,s)] ds,a,b)+ int(Cc2(s)[x_2ss(t,s)] ds,a,b) + [D21]w(t) + [D22]u(t)
%                       x_2(b)             [x_2(t,s)]                     [x_3s(t,s)]
%                       x_3(a)             [x_3(t,s)]
%                       x_3(b)
%                       x_3s(a)
%                       x_3s(b)]
%
%
% % Boundary conditions are of the form
% % B[x_2(a)
% %   x_2(b)
% %   x_3(a)
% %   x_3(b)
% %   x_3s(a)
% %   x_3s(b)]= Bx xo(t)+ Bw w(t)+ Bu u(t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-Defined Inputs for specifying the dynamics:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: All dimensions (nx,nw,nu,nz,ny) are assumed to be 0 by default.
% Specifically, (all optional)
% nx     -  number of ODE states
% nw     -  number of disturbances
% nu     -  number of controlled inputs
% nz     -  number of regulated outputs
% ny     -  number of observed outputs
%
% NOTE: Any matrix with a 0-dimension should be ommitted
%
% PDE Terms (required)
% n0     -  number of undifferentiated PDE states (Required)
% n1     -  number of single differentiated PDE states (Required)
% n2     -  number of double differentiated PDE states (Required)
% a,b    -  \R x \R - interval of the domain of the spatial variable - s \in [a,b] (Required)
% B      -  matrix of dimension n1+2*n2 x n0+n1+n2 (Required, must have row rank n1+2n2 and contain no prohibited boundary conditions)
% A0(s)  -  matrix function of s of dimension n0+n1+n2 x n0+n1+n2
% A1(s)  -  matrix function of s of dimension n0+n1+n2 x n1+n2
% A2(s)  -  matrix function of s of dimension n0+n1+n2 x n2
%%%%%% Input coupling (optional)
% Bw     -  matrix of dimension n1+2*n2 x nw (must be full row rank)       - effect of disturbance on Boundary Conditions
% Bu     -  matrix of dimension n1+2*n2 x nu (must be full row rank)       - effect of the input on the boundary conditions
% B21(s) -  polynomial matrix in pvar s of dimension (n0+n1+n2) x nw       - Distributed effect of disturbance in the domain of the PDE
% B22(s) -  polynomial matrix in pvar s of dimension (n0+n1+n2) x nu       - Distributed effect of input in the domain of the PDE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Dynamics (all optional)
% no     -  number of ODE states
% A      - matrix of dimension nx x nx                                     - ODE dynamics
% B11    -  matrix of dimension nx x nw                                    - effect of disturbance on ODE state
% B12    -  matrix of dimension nx x nu                                    - effect of input on ODE state
%%%%%%%% Coupling of ODE to PDE
% Bx    -  matrix of dimension n1+2*n2 x nx                                - Effect of ODE state on boundary conditions
% E(s)   - pvar matrix in variable s of dimension n0+n1+n2 x nx            - Distributed effect of ODE in the domain of the PDE
%%%%%%%% Coupling of PDE to ODE
% E0     -  matrix of dimension nx x 2*n1+4*n2                             - Effect of boundary values of distributed states on ODE states
% Ea     -  polynomial matrix in pvar s of dimension nx x (n0+n1+n2)       - kernel of integration for effect of distributed states on ODE states
% Eb     -  polynomial matrix in pvar s of dimension nx x n1+n2            - kernel of integration for effect of 1st-order derivatives of distributed states on ODE states
% Ec     -  polynomial matrix in pvar s of dimension nx x n2               - kernel of integration for effect of 2nd-order derivatives of distributed states on ODE states
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disturbance-related inputs (all optional)
% D21 (see observed outputs), D11 (see regulated output)
% Bw,B21 (See PDE dynamics); B11 (See ODE dynamics)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control-related inputs (all optional)
% D22 (see observed outputs); D12 (see regulated output)
% Bu,B22 (See PDE dynamics); B12 (See ODE dynamics)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regulated outputs (all optional)
% nz     -  number of regulated outputs
% C1     -  matrix of dimension nz x nx                                    - Effect of ODE states on regulated output
% D11    -  matrix of dimension nz x nw                                    - effect of disturbance on regulated output (avoid for H2-optimal control problems)
% D12    -  matrix of dimension nz x nu                                    - effect of control input on regulated output (Recommended for realistic controllers)
% C10    -  matrix of dimension nz x 2*n1+4*n2                             - Effect of boundary values of distributed states on regulated output
% Ca1(s) -  polynomial matrix in pvar s of dimension nz x (n0+n1+n2)       - kernel of integration for effect of distributed states on regulated outputs
% Cb1(s) -  polynomial matrix in pvar s of dimension nz x n1+n2            - kernel of integration for effect of 1st-order derivatives of distributed states on regulated outputs
% Cc1(s) -  polynomial matrix in pvar s of dimension nz x n2               - kernel of integration for effect of 2nd-order derivatives of distributed states on regulated outputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observed outputs (all optional)
% ny     -  number of observed outputs
% C2     -  matrix of dimension ny x nx                                    - Effect of ODE states on observed output
% D21    -  matrix of dimension ny x nw                                    - effect of disturbance on observed output (e.g. sensor noise)
% D22    -  matrix of dimension ny x nu                                    - effect of control input on observed output (rare)
% C20    -  matrix of dimension ny x 2*n1+4*n2                             - Effect of boundary values of distributed states on observed outputs
% Ca2(s) -  polynomial matrix in pvar s of dimension ny x (n0+n1+n2)       - kernel of integration for effect of distributed states on observed outputs
% Cb2(s) -  polynomial matrix in pvar s of dimension ny x n1+n2            - kernel of integration for effect of 1st-order derivatives of distributed states on observed outputs
% Cc2(s) -  polynomial matrix in pvar s of dimension ny x n2               - kernel of integration for effect of 2nd-order derivatives of distributed states on observed outputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Library Call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To test PIETOOLS on one of the standard example problems, consult the
% library file examples_PDE_library_PIETOOLS.m and uncomment the example
% you wish to test. Otherwise, comment the following line and enter your
% own PDE.

examples_PDE_library_PIETOOLS

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

% NOTE: At present, PIETOOLS does not support inputs at the boundary for
% solving the Hinf optimal control problem. Support for this option will be
% included in a future release. Use of the Hinf_control=1 toggle in this
% case will generate an error.


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

sosineq_on=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strictness of Inequalities:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eppos = 1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
eppos2 = 1*1e-6;   % Positivity of Lyapunov Function with respect to spatially distributed states
epneg = 0*1e-5;    % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Use of heavy settings may be recommended for computing Hinf gain or
% optimal controllers/estimators
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% settings_PIETOOLS_veryheavy
settings_PIETOOLS_heavy
%  settings_PIETOOLS_light
% settings_PIETOOLS_stripped
% settings_PIETOOLS_extreme
% settings_PIETOOLS_custom


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END: Linear PI Inequality (LPI) solver settings:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Convert your ODE-PDE system to a PIE
convert_PIETOOLS_PDE;

%%% Construct and solve the LPI corresponding to the problem you specified earlier
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

% NOTE: if multiple executives are called, prog only contains data from the
% most recent executive!
