% DEMO2_volterra_operator_norm.m
% See Chapter 11.2 of the manual for a description.
%
% This document illustrates how an upper bound on the norm of the Volterra 
% integral operator can be computed in PIETOOLS.
% Volterra integral operator
%   (Top*x)(s) = int(x(r),dr,a,s),        s in [a,b]
% Optimization Problem (LPI)
% min   γ,
% s.t   T'*T ≤ γ
%
% This example is also included in the paper (page 5, Demonstration 2)
% link: https://arxiv.org/pdf/1910.01338.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO2
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
% DJ, 11/19/2024: Simplify demo (remove lines of code where possible);

clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare the operators of interest

% % Define Volterra operator
% %  (Top*x)(s) = int_{a}^{s} x(r) dr    s in [a,b]
a=0;    b=1; 
opvar Top;
Top.R.R1 = 1; Top.I = [a,b];                   

% =============================================
% === Declare the LPI

% % Initialize LPI program
prob = lpiprogram(Top.vars,Top.I);

% % Declare decision variables:
% %   gam \in \R
[prob,gam] = lpidecvar(prob,'gam');

% % Set inequality constraints
% %   Top'*Top-gam <= 0
opts.psatz = 1;                             % Allow Top'*Top-gam>0 outside of [a,b]
prob = lpi_ineq(prob,gam-Top'*Top,opts); % lpi_ineq(prob,Q) enforces Q>=0

% % Set objective function:
% %   min gam
prob = lpisetobj(prob, gam);

% % Solve LPI and retrieve solution
prob = lpisolve(prob);
operator_norm = sqrt(double(lpigetsol(prob,gam)));


echo off

fprintf(['\n If successful, ',num2str(operator_norm),' is an upper bound on the norm of the Volterra integral operator.\n'])
fprintf([' The exact operator norm of the Volterra integral operator is 2/pi=',num2str(2/pi),'.\n']);