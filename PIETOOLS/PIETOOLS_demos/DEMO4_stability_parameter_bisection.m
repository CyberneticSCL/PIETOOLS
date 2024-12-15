% DEMO4_stability_parameter_bisection.m
% See Chapter 11.4 of the manual for a description.
%
% This document illustrates how PIETOOLS can be used to find the maximal
% value of a parameter for which a PDE is stable, using bisection.
% Specifically, we test stability of the reaction-diffusion PDE
%   x_{t}(t,s) = lam*x(t,s) + x_{ss}(t,s);
%       x(t,0) = x(t,1) = 0;
% The PDE is stable when lam <= pi^2 = 9.8696 (Ahmadi 2015).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO4
%
% Copyright (C)2024  PIETOOLS Team
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
% MP, SS, DJ, 2022: Initial coding;
% DJ, 10/20/2024: Update to use new LPI programming functions;
% DJ, 11/19/2024: Simplify demo (remove lines of code where possible, use 
%                   'sys' structure to declare PDE);

clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare parameters

% Declare independent variables (time and space)
pvar t s

% Set bisection limits for lam
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