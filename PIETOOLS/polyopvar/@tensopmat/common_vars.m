function [A_out,B_out,old2new_idcs,new2old_idcs] = common_vars(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_VARS(A_IN,B_IN) takes two tensopmat objects
% and expresses them in terms of the same spatial variables
%
% INPUTS
% - A_in, B_in:     'tensopvar' objects
%
% OUTPUTS
% - A_out, B_out:   'tensopmat' objects representing the same tensor-PI
%                   operators as A_in, B_in, respectively, but now 
%                   expressed in terms of the same variables, so that
%                       A_out.vars = B_out.vars;    A_out.dom = B_out.dom;
% - old2new_idcs:   N x 1 array of indices specifying how the new variables
%                   can be expressed in terms of the concatenated list of
%                   variables:
%                       A_out.vars = vars(old2new_idcs,:)
%                   where vars = [A_in.vars; B_in.vars];
% - new2old_idcs:   N1+N2 x 1 array or indices specifying how the old
%                   variables can be recovered from the new variables:
%                       [A_in.vars; B_in.vars] = A_out.vars(new2old_idcs,:)
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
% DJ, 04/22/2026: Initial coding

% Extract the variables appearing in each tensor-PI operator
vars_A = A_in.vars;         vars_B = B_in.vars;
dom_A = A_in.dom;           dom_B = B_in.dom;
depmat1_A = A_in.depmat1;   depmat1_B = B_in.depmat1;
depmat2_A = A_in.depmat2;   depmat2_B = B_in.depmat2;

% Deal with case of constant coefficient
if all(depmat1_A==0) && all(depmat2_A==0)
    vars_A = vars_B;
    dom_A = dom_B;
elseif all(depmat1_B==0) && all(depmat2_B==0)
    vars_B = vars_A;
    dom_B = dom_A;
end

% Combine the variable names into a list of unique variables
varname_A = pvar2varname(vars_A(:,1));    nvar_A = numel(varname_A);
varname_B = pvar2varname(vars_B(:,1));
if ~isempty(A_in.ops) && isa(A_in.ops{1},'intop') && ~isequal(varname_A,varname_B)
    error("The variables defining the different 'intop' objects do not match, this may lead to errors.")
end
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
depmat1_A_new = zeros(size(depmat1_A,1),nvar);
depmat1_A_new(:,var_idcs_A) = depmat1_A;
depmat2_A_new = zeros(size(depmat2_A,1),nvar);
depmat2_A_new(:,var_idcs_A) = depmat2_A;
depmat1_B_new = zeros(size(depmat1_B,1),nvar);
depmat1_B_new(:,var_idcs_B) = depmat1_B;
depmat2_B_new = zeros(size(depmat2_B,1),nvar);
depmat2_B_new(:,var_idcs_B) = depmat2_B;

% Declare the tensor-PI operators in terms of the same variables
A_out = A_in;                   B_out = B_in;
A_out.vars = vars_AB;           B_out.vars = vars_AB;
A_out.dom = dom_AB;             B_out.dom = dom_AB;
A_out.depmat1 = depmat1_A_new;  B_out.depmat1 = depmat1_B_new;
A_out.depmat2 = depmat2_A_new;  B_out.depmat2 = depmat2_B_new;

end