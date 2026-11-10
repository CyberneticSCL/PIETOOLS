function [A_out,B_out] = common_vars(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_VARS(A_IN,B_IN) takes two 'tensopvar' objects and
% expresses them in terms of the same spatial variables.
%
% INPUTS
% - A_in, B_in:     'tensopvar' objects;
%
% OUTPUTS
% - A_out, B_out:   'tensopvar' objects representing the same tensor-PI
%                   operators as A_in, B_in, respectively, but expressed
%                   in terms of the same variables, so that
%                       A_out.vars = B_out.vars;    A_out.dom = B_out.dom;
%
% NOTES
% - The dependency arrays depmat.in/.out of A_out, B_out are extended to
%   the combined list of variables.
% - Variables shared by A_in and B_in must have the same domain.
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
% DJ, 10/15/2026: Initial coding;

% Extract the variables appearing in each tensor-PI operator
vars_A = A_in.vars;             vars_B = B_in.vars;
dom_A = A_in.dom;               dom_B = B_in.dom;
depmat_A = A_in.depmat;         depmat_B = B_in.depmat;

% Establish unique input and output variables
[vars1,dom1,depmat1_A_new,depmat1_B_new] = combine_vars(vars_A.out,dom_A.out,depmat_A.out,vars_B.out,dom_B.out,depmat_B.out);
[vars2,dom2,depmat2_A_new,depmat2_B_new] = combine_vars(vars_A.in,dom_A.in,depmat_A.in,vars_B.in,dom_B.in,depmat_B.in);

% Declare the tensor-PI operators in terms of the same variables
A_out = A_in;                       B_out = B_in;
A_out.vars.out = vars1;             B_out.vars.out = vars1;
A_out.vars.in = vars2;              B_out.vars.in = vars2;
A_out.dom.out = dom1;               B_out.dom.out = dom1;
A_out.dom.in = dom2;                B_out.dom.in = dom2;
A_out.depmat.out = depmat1_A_new;   B_out.depmat.out = depmat1_B_new;
A_out.depmat.in = depmat2_A_new;    B_out.depmat.in = depmat2_B_new;

end



function [vars,dom,depmatA_new,depmatB_new] = combine_vars(varsA,domA,depmatA,varsB,domB,depmatB)
% Establish a unique list of variables and associated domains given the
% variables [varsA;varsB] and domains [domA;domB]. Update the dependency
% arrays to match this full list of variables.

% Combine the variable names into a list of unique variables
nvarA = numel(varsA);
if isequal(varsA,varsB)
    varsAB = varsA;
    old2new_idcs = (1:nvarA)';
    new2old_idcs = [1:nvarA,1:nvarA]';
else
    [varsAB,old2new_idcs,new2old_idcs] = unique([varsA,varsB]);
end
var_idcs_A = new2old_idcs(1:nvarA);
var_idcs_B = new2old_idcs(nvarA+1:end);
nvars = numel(varsAB);

% Establish unique variables and associated domains
vars_full = [varsA, varsB];
vars = vars_full(old2new_idcs);
dom_full = [domA; domB];
dom = dom_full(old2new_idcs,:);
if any(any(dom(new2old_idcs,:)~=dom_full))
    error("Identical variables in the operators must have the same domain.")
end

% Augment depmats to account for merging of variables
depmatA_new = zeros(size(depmatA,1),nvars);
depmatA_new(:,var_idcs_A) = depmatA;
depmatB_new = zeros(size(depmatB,1),nvars);
depmatB_new(:,var_idcs_B) = depmatB;

end