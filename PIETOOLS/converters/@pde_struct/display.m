function display(PDE,name)
% display(PDE,name) 
% Displays in the command window the PDE defined by the
% pde_struct object "PDE", with name "name".
%
% INPUTS:
% PDE:  A pde_struct class object, defining a PDE in a manner as
%       outlined in "initialize_PIETOOLS_PDE".
% name: str, specifying the name of the PDE. (not currently used...)
% 
% OUTPUTS:
% Command line display of PDE, as well as the output equations and BCs.
% 
% NOTES:
% - Although spatial variable names may be specified in the PDE structure,
%   these will not be used in the display. Instead, s_i will be used to
%   denote primary variables, and phi_i to denote dummy variables.
% - Coefficient matrices term{j}.C will not be shown, unless they are a
%   constant scalar value. If they are matrix-valued and/or polynomial, "C"
%   will be used to denote the factor, to avoid the expression from getting
%   too messy.
% - Display is only possible for initialized PDEs, as the fields x_tab
%   through BC_tab constructed in "initialize_PIETOOLS_PDE" are used to
%   establish variable dependence of the different states, inputs, and
%   outputs. An unfinished PDE will likely produce issues in the
%   initialization function, so only complete systems can be displayed.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07/08/2022

if nargin<2
    name = '';
end

try PDE = initialize_PIETOOLS_PDE(PDE,true);
    display_PDE(PDE,name);
catch me
    disp(PDE);
end

return