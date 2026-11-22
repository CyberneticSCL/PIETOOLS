%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_FI.m     PIETOOLS 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Started: This file is meant to be run from the Matlab editor
% window. It computes the constant of a functional inequality from the
% library examples_FI_library_PIETOOLS.
%
% Step 1: Choose an inequality and its options
% Step 2: Compute its constant as the solution of an LPI
%
% Each inequality of the library holds on [0,b]. It is written as an
% operator inequality in the free state y = x_s (and, where needed, the
% value x(0)), so its constant C solves the LPI
%       min C   s.t.   (quadratic form of operators with C) >= 0.
% The library holds 17 inequalities: Poincare-type (1-5), trace (6-9),
% Hardy (10), operator norms (11-14) and eigenvalue bounds (15-17); the
% table of contents is in its help (help examples_FI_library_PIETOOLS).
%
% INPUT
% - none; Step 1 selects example 1, the Wirtinger inequality
%   int x^2 <= C int x_s^2 with x(0) = x(b) = 0, on [0,1].
%
% OUTPUT
% - none; the workspace receives C, the constant the LPI returns, and info
%   (fields name, constant, sharp, paper, extra), as
%   examples_FI_library_PIETOOLS returns them. The library prints C next
%   to the sharp constant where one is known.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 PIETOOLS Team
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding MMP, 11/22/2026: the entry script of the functional-
%   inequality library, as PIETOOLS_PDE, PIETOOLS_DDE and PIETOOLS_DDF are
%   of the other three libraries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Choose an inequality and its options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The index of the inequality, 1 through 17 (help examples_FI_library_PIETOOLS)
example = 1;
% Options, each a 'name=value' string: b (the interval [0,b], default 1);
% dD, dw (degrees of the diagonal and off-diagonal parts of the
% certificate); t (the point of example 8); eps (the weight of example
% 10); tmin (the range of t in example 9); modes (eigenvalues of example
% 15, eigenfunctions of example 17); plot (1 plots the eigenfunctions of
% examples 15 and 17, 0 does not)
options = {'b=1'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Compute the constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[C,info] = examples_FI_library_PIETOOLS(example,options{:});
