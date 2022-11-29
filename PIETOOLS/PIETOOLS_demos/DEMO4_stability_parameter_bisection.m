% DEMO4_stability_parameter_bisection.m
% See Chapter 11.4 of the manual for a description.
%
% This document illustrates how PIETOOLS can be used to find the maximal
% value of a parameter for which a PDE is stable, using bisection.

% Reaction-diffusion PDE
% \dot{x}(t,s) = lam*x(t,s) + x_{ss}(t,s);
% x(t,0) = x(t,1) = 0;
% Stable when lam <= pi^2 = 9.8696 (Ahmadi 2015).

%%
clc; clear;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

%%% Set bisection limits for lam.
lam_min = 0;        lam_max = 20;
lam = 0.5*(lam_min + lam_max);
n_iters = 8;

%%% Initialize a PDE structure.
a = 0;  b = 1;
pvar s
pde_struct PDE;
PDE.x{1}.vars = s;
PDE.x{1}.dom = [a,b];

% Set the PDE \dot{x}(t,s) = lam*x(t,s) + x_{ss}(t,s);
PDE.x{1}.term{1}.C = lam;
PDE.x{1}.term{2}.D = 2;

% Set the BCs x(t,a) = x(t,b) = 0;
PDE.BC{1}.term{1}.loc = a;
PDE.BC{2}.term{1}.loc = b;

%%% Initialize settings for solving the LPI
settings = lpisettings('veryheavy');
if strcmp(settings.sos_opts.solver,'sedumi')
    settings.sos_opts.params.fid = 0;   % Suppress output in command window
end

%%% Perform bisection on the value of lam
for iter = 1:n_iters
    % Update the value of lam in the PDE.
    PDE.x{1}.term{1}.C = lam;
    
    % Update the PIE.
    PIE = convert(PDE,'pie');
    
    % Re-run the stability test.
    fprintf(['\n',' --- Running the stability test with lam = ',num2str(lam),' ---\n'])
    prog = PIETOOLS_stability(PIE,settings);

    % Check if the system is stable
    if prog.solinfo.info.dinf || prog.solinfo.info.pinf || abs(prog.solinfo.info.feasratio - 1)>0.3
        % Stability cannot be verified, decreasing the value of lam...
        lam_max = lam;
        lam = 0.5*(lam_min + lam_max);
    else
        % The system is stable, trying a larger value of lam...
        lam_min = lam;
        lam = 0.5*(lam_min + lam_max);
    end    
end
%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off

fprintf(['\n Stability of the system could be verified for lam<=',num2str(lam_min),'.\n'])
fprintf([' An analytic bound on lam guaranteeing stability is pi^2=',num2str(pi^2),'.\n']);

% @article{valmorbida2015stability,
%   title={Stability analysis for a class of partial differential equations via semidefinite programming},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   journal={IEEE Transactions on Automatic Control},
%   volume={61},
%   number={6},
%   pages={1649--1654},
%   year={2015},
%   publisher={IEEE}
% }