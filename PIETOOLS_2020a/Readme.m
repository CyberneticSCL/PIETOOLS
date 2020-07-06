% List of PIETOOLS functions and descriptions
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - README
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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



% Directory Structure of PIETOOLS:
% - converters: scripts for initialization and conversion of DDE,DDF/PDE to PIE
% - executives: scripts for solving LPIs for stability, hinf gain, observer and controller design
% - multipoly: multipolynomial toolbox
% - opvar: functions and methods related to opvar class 
% - settings: standard settings for LPI optimization problems
% - SOSTOOLS.303: SOSTOOLS toolbox
% - SOSTOOLS_ADD: modified files of SOSTOOLS
% - GetStarted_DOCS_DEMOS: demonstration code and document


% Files in PIETOOLS:
% - converters (directory)
% -- convert_PIETOOLS_DDE: checks DDE definition is accurate and converts to PIE
% -- convert_PIETOOLS_DDF: checks DDF definition is accurate and converts to PIE
% -- convert_PIETOOLS_PDE: checks PDE definition is accurate and converts to PIE
% -- initialize_PIETOOLS_DDE: checks DDE definition is accurate
% -- initialize_PIETOOLS_DDF: checks DDF definition is accurate
% -- initialize_PIETOOLS_PDE: checks PDE definition is accurate

% - executives (directory)
% -- executive_PIETOOLS_Hinf_control: Initializes and solves hinf-optimal controller LPI for given PIE
% -- executive_PIETOOLS_Hinf_estimator: Initializes and solves hinf-optimal estimator LPI for given PIE
% -- executive_PIETOOLS_Hinf_gain: Initializes and solves hinf-gain LPI for given PIE
% -- executive_PIETOOLS_Hinf_gain_dual: Initializes and solves hinf-gain dual LPI for given PIE
% -- executive_PIETOOLS_stability: Initializes and solves stability LPI for given PIE
% -- executive_PIETOOLS_stability_dual: Initializes and solves stability dual LPI for given PIE

% - GetStarted_DOCS_DEMOS (directory)
% -- GetStarted_DEMO: An example of LPI typically solved using PIETOOLS
% -- PIETOOLS_2020a_How_to_get_started: A pdf explaining the steps in solving LPIs
% -- poincare_inequiality_DEMO: Solves for poincare constant in 1D
% -- conversion_pde2pie_DEMO: Initializes pde and converts to pie
% -- volterra_operator_norm_DEMO: Solves for operator norm of volterra integral operator

% - opvar (directory)
% -- @opvar (directory)
% ------ ctranspose: Finds adjoint of opvar object
% ------ degbalance: Finds degrees of poslpivar needed to solve an LPI
% ------ display: Overrides default display format of opvar class
% ------ eq: Tests for equality of opvar objects; Returns logical value
% ------ getdeg: Finds least and highest polynomial degrees in opvar elements
% ------ horzcat: Performs horizontal concatenation of opvar objects
% ------ isempty_opvar: Tests if an opvar object is empty; Returns logical value
% ------ isvalid: Tests if an opvar object has consistent dimensions and standard independent variables
% ------ minus: Subtracts two opvar objects
% ------ mtimes: Composes two opvar objects
% ------ plus: Adds two opvar objects
% ------ op_slice: Uses index method to generate slices of opvar components
% ------ uminus: Unitary minus operator for opvar object
% ------ uplus: Unitary plus operator for opvar object
% ------ vertcat: Performs vertical concatenation of opvar objects
% -- getsol_lpivar: Finds the solution of lpivar from the solved LPI problem
% -- inv_opvar: Numerical inverse of an invertible opvar object (prototype, use with caution)
% -- lpi_eq: Adds PI-valued equality constraints to LPI problem
% -- lpi_ineq: Adds PI-valued inequality constraints to LPI problem
% -- lpivar: Initializes an opvar variable to be used in LPI
% -- monomialdiff: Finds the set difference in monomial set used to build two polynomials
% -- opvar_postest: Numerically tests if an opvar object is sign definite
% -- poslpivar: Initializes a positive opvar variable to be used in LPI
% -- rand_opvar: Initializes a random opvar object
% -- subs_op: Substitutes a variable with another specified value
% -- show: Alternative display function for opvar objects that prints selected components

% - settings (directory)
% -- settings_PIETOOLS_custom: Settings with customizable optimization parameters
% -- settings_PIETOOLS_extreme: Settings that produce very small-size optimization problem, low degree polynomials, highly sparse lpivars
% -- settings_PIETOOLS_heavy: Settings that produce large-size optimization problem, high degree polynomials, non sparse lpivars
% -- settings_PIETOOLS_light: Settings that produce small-size optimization parameters, low degree polynomials, moderately-sparse lpivars
% -- settings_PIETOOLS_veryheavy: Settings that produce very large-size optimization problem, very high degree polynomials, non sparse lpivars
% -- settings_PIETOOLS_stripped: Settings that produce smallest-size optimization parameters, low degree polynomials, highly sparse lpivars

% - (root directory)
% -- examples_dde_library_PIETOOLS: Script containing examples of DDEs
% -- examples_ddf_library_PIETOOLS: Script containing examples of DDFs
% -- examples_pde_library_PIETOOLS: Script containing examples of PDEs
% -- PIETOOLS_DDE: solver that initializes DDE, converts to PIE and solves a chosen LPI
% -- PIETOOLS_DDF: solver that initializes DDF, converts to PIE and solves a chosen LPI 
% -- PIETOOLS_PDE: solver that initializes PDE, converts to PIE and solves a chosen LPI
% -- Readme: The file you are currently reading!! 

