function [A_out,B_out] = common_vars(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_VARS(A_IN,B_IN) takes two tensopvar objects
% and expresses them in terms of the same spatial variables
%
% INPUTS
% - A_in, B_in:     'tensopvar' objects
%
% OUTPUTS
% - A_out, B_out:   'tensopvar' objects representing the same tensor-PI
%                   operators as A_in, B_in, respectively, but now 
%                   expressed in terms of the same variables, so that
%                       A_out.var1 = B_out.var1;    A_out.dom1 = B_out.dom1;
%                       A_out.var2 = B_out.var2;    A_out.dom2 = B_out.dom2;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - common_vars
%
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
% DJ, 04/15/2026: Initial coding

% Extract the variables appearing in each tensor-PI operator
vars_A = A_in.vars;         vars_B = B_in.vars;
dom_A = A_in.dom;           dom_B = B_in.dom;
depmat1_A = A_in.depmat1;   depmat1_B = B_in.depmat1;
depmat2_A = A_in.depmat2;   depmat2_B = B_in.depmat2;

% Combine the variable names into a list of unique variables
varname_A = pvar2varname(vars_A(:,1));    nvar_A = numel(varname_A);
varname_B = pvar2varname(vars_B(:,1));
[varname_AB,old2new_idcs,new2old_idcs] = unique([varname_A;varname_B]);
var_idcs_A = new2old_idcs(1:nvar_A);
var_idcs_B = new2old_idcs(nvar_A+1:end);
nvar = numel(varname_AB);

% Establish unique variables and associated domains
vars_AB = [vars_A; vars_B];
vars_AB = vars_AB(old2new_idcs,:);
dom_AB = [dom_A; dom_B];
dom_AB = dom_AB(old2new_idcs,:);

% Augment depmats to account for merging of variables
depmat1_A_new = zeros(nvar,size(depmat1_A,2));
depmat1_A_new(var_idcs_A,:) = depmat1_A;
depmat2_A_new = zeros(nvar,size(depmat2_A,2));
depmat2_A_new(var_idcs_A,:) = depmat2_A;
depmat1_B_new = zeros(nvar,size(depmat1_B,2));
depmat1_B_new(var_idcs_B,:) = depmat1_B;
depmat2_B_new = zeros(nvar,size(depmat2_B,2));
depmat2_B_new(var_idcs_B,:) = depmat2_B;

% Declare the tensor-PI operators in terms of the same variables
A_out = A_in;                   B_out = B_in;
A_out.vars = vars_AB;           B_out.vars = vars_AB;
A_out.dom = dom_AB;             B_out.dom = dom_AB;
A_out.depmat1 = depmat1_A_new;  B_out.depmat1 = depmat1_B_new;
A_out.depmat2 = depmat2_A_new;  B_out.depmat2 = depmat2_B_new;

end