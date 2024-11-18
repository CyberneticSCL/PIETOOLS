% DEMO4_stability_parameter_bisection.m
% See Chapter 11.4 of the manual for a description.
%
% This document illustrates how PIETOOLS can be used to find the maximal
% value of a parameter for which a PDE is stable, using bisection.
% Specifically, we test stability of the reaction-diffusion PDE
%   x_{t}(t,s) = lam*x(t,s) + x_{ss}(t,s);
%       x(t,0) = x(t,1) = 0;
% The PDE is stable when lam <= pi^2 = 9.8696 (Ahmadi 2015).

clc; clear; clear stateNameGenerator;
echo on
%% %%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%
pvar s t

% =============================================
% === Declare parameters

% Set bisection limits for lam.
lam_min = 0;        lam_max = 20;
lam = 0.5*(lam_min + lam_max);
n_iters = 8;

% Initialize settings for solving the LPI
settings = lpisettings('heavy');
if strcmp(settings.sos_opts.solver,'sedumi')
    settings.sos_opts.params.fid = 0;   % Suppress output in command window
end

%%% Perform bisection on the value of lam
for iter = 1:n_iters
    % =============================================
    % === Declare the operators of interest

    % Declare a PDE 
    PDE = sys();
    x = state('pde');
    PDE = addequation(PDE, [diff(x,t)==diff(x,s,2)+lam*x;
                            subs(x,s,0)==0;
                            subs(x,s,1)==0]);

    % Convert to PIE.
    PIE = convert(PDE,'pie');

    % =============================================
    % === Declare the LPI

    % Run pre-defined stability executive
    fprintf(['\n',' --- Running the stability test with lam = ',num2str(lam),' ---\n'])
    prog = lpiscript(PIE,'stability',settings);

    % Check if the system is stable
    if prog.solinfo.info.dinf || prog.solinfo.info.pinf || abs(prog.solinfo.info.feasratio - 1)>0.3
        % Stability cannot be verified --> decrease value of lam...
        lam_max = lam;
        lam = 0.5*(lam_min + lam_max);
    else
        % The system is stable --> try larger value of lam...
        lam_min = lam;
        lam = 0.5*(lam_min + lam_max);
    end    
end
%% %%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
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