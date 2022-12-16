% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - README
%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% This file is packaged with PIETOOLS and contains License information, a
% list of all the files packaged with PIETOOLS and a brief description of
% the functions and scripts. If any of the following files are missing from
% the PIETOOLS folder, the toolbox may not function as expected. In that
% case, reinstall using the install script or contact sshivak@asu.edu or
% mpeet@asu.edu for support.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What is PIETOOLS?
% PIETOOLS is a free MATLAB toolbox for formulating and solving Linear PI 
% Inequalities (LPIs) programs. PIETOOLS can be used to define 1D and 2D PI 
% operators, declare 1D and 2D PI operator variables (postive semidefinite 
% or indefinite), add operator inequality constraints, and solve LPI 
% optimization problems. The interface is inspired by YALMIP and the program 
% structure is based on that used by SOSTOOLS. By default the LPIs are solved
% using SeDuMi.
%
% What is a PIE?
% PIE stands for Partial Integral Equation and is an alternative 
% representation for many commonly encountered classes of systems, including 
% Ordinary Differential Equations (ODEs), Partial-Differential Equations 
% (PDEs), Delay Differential Equations (DDEs), and Differential-Difference 
% Equations (DDFs).
% 
% The cool thing about PIEs is that, unlike PDEs and DDEs (which have lame 
% boundary conditions, unbounded operators, and continuity constraints), PIEs
% are defined by the very slick linear algebra of 3/4 PI operators. This 
% feature makes PIEs the representation of choice if you want to do anything 
% computational with your beam equation, network model, reaction-diffusion 
% equations, etc.
% 
% Now, you may be wondering if you are going to lose anything by switching 
% your PDE/DDE/DDF to a PIE. No! That would be awful. You may have been hurt 
% in the past by people wrecking your lovely PDE/DDE/DDF using such barbaric 
% tools as approximation via discretization, projection, mollification, 
% regularization or Pade. However, let us assure you that using PIEs is 
% completely safe. The PIE representation of a PDE/DDE/DDF is exact. The 
% solutions are one-to-one, only the tools used for representation have 
% improved.
% 
% What can I do with PIETOOLS?
% The use of PIETOOLS is organized into 4 categories. These are:
% 1) Manipulation of 1D and 2D PI operators.
% 2) Constructing and Solving Linear PI Inequalities (LPIs).
% 3) Converting PDEs/DDEs/DDFs to PIEs.
% 4) Solving LPIs for Analysis and Control of PIEs.
% 5) Numerically simulating PDE, DDE, and PIE solutions.