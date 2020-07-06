% This document illiustrates the use of PIETOOLS to verify the 
% stability of a reaction-diffusion equation.

% You may like to consult PIETOOLS_2020a__How_to_get_started.pdf
% to know more about each line of the code. 

% PDE Model
% PDE:     \dot(x)(t,s) = \lambda x(t,s) +  x_ss(t,s),   a<s<b
% BCs:     x(t,a) =0,  x(t,b) =0
% Domain:  a=0, b=1

clc; clear;
%%
%%% STEP 1: Specify the model
% Number of state variables based on their derivatives
n0=0;n1=0;n2=1; 

% Define PDE parameters
lambda = 10; %10
A0 = lambda; A1 = 0; A2=1;

% Define the Domain [a,b]
a=0;b=1;

% Define Boundary conditions in the form B[x(a); x(b); x_s(a) x_s(b)]= 0
B=[1 0 0 0;
	0 1 0 0];




%%
%%% STEP 2: Construct and View the PIE Representation

% Define the polynomial variable s, theta. DO NOT CHANGE THEM!
pvar s theta;  

% Use conversion script to convert PDE to PIE
convert_PIETOOLS_PDE;

% View the corresponding PI Operators
Aop

Top

%%
%%% STEP 3: Solve The Stability Problem

% Initialize the optimization program with independepent variables 
prog = sosprogram([s; theta]);

% Declare an unknown positive PI variable Pop
[prog, Pop] = poslpivar(prog, [0 ,1],[0 1]);


Pop = Pop+.0001; 

% Declare the LPI corresponding to stability test

options.psatz = 1;
prog = lpi_ineq(prog,-[Top'*Pop*Aop+Aop'*Pop*Top], options);

% Solve the feasibility test
prog = sossolve(prog);

%%
%%% STEP 4: Display the result

% Run diagnostics for feasibility and numerical error
if norm(prog.solinfo.info.feasratio-1)<=.1 
    disp('The System of equations was successfully solved.')
else
    disp('The System is of equations not solved.')
end


