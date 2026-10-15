function [A_out,B_out] = common_intmat(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_INTMAT(A_IN,B_IN) expresses the 'intop' objects
% A_IN and B_IN in terms of the same set of integrals.
%
% INPUTS
% - A_in:   'intop' object;
% - B_in:   'intop' object of the same dimensions as A_in, acting on the
%           same function space, so that A_in.vars = B_in.vars,
%           A_in.depmat = B_in.depmat, and A_in.dom = B_in.dom;
%
% OUTPUTS
% - A_out:  'intop' object representing the same functional as A_in,
%           expressed in terms of the same integrals as B_out, so that
%           A_out.intmat = B_out.intmat;
% - B_out:  'intop' object representing the same functional as B_in,
%           expressed in terms of the same integrals as A_out;
%
% NOTES
% - A_in and B_in are first expressed in terms of the same monomials and
%   decision variables, if they are not yet.
% - An integral that appears in only one of the functionals is added to the
%   other with zero coefficients.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - common_intmat
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

% Check that both objects are properly declared.
if (~A_in.chkval && ~isvalid(A_in)) || (~B_in.chkval && ~isvalid(B_in))
    error("One of the input 'intop' objects is not properly declared.")
end

% Check that the functionals act on the same space
if any(size(A_in)~=size(B_in))
    error("Functionals must be of the same matrix dimensions.")
elseif ~isequal(A_in.vars,B_in.vars)
    error("Functionals must be defined in terms of the same variables.")
elseif ~isequal(A_in.dom,B_in.dom)
    error("Functionals must be defined on the same spatial domain.")
elseif ~isequal(A_in.depmat,B_in.depmat)
    error("Functionals must act on the same function space.")
end

% Express the functionals in terms of the same monomials and decision variables
[A_in,B_in] = sync_basis(A_in,B_in);

% Make sure intmats for both functionals do not contain duplicates
A_in = unique_intmat(A_in);
B_in = unique_intmat(B_in);

% Extract the intervals of integration and associated kernels
intmat1 = A_in.intmat;  params1 = A_in.params;
intmat2 = B_in.intmat;  params2 = B_in.params;

% Combine intmats into unique matrices
N = numel(intmat1);
if N==0 || isequal(intmat1,intmat2)
    % No variables, or the same integrals: already in terms of the same
    % integrals, with the parameters in the same order
    A_out = A_in;   B_out = B_in;
    return
end
intmat_new = cell(1,N);     nterms_arr = zeros(1,N);
ridcs1_arr = cell(1,N);     ridcs2_arr = cell(1,N);
cidcs1_arr = cell(1,N);     cidcs2_arr = cell(1,N);
for i=1:numel(intmat1)      % could be parallelized
    % Get full list of intervals for variable i
    intmat1_i = intmat1{i};     nr1 = size(intmat1_i,1);
    intmat2_i = intmat2{i};
    intmat_full = [intmat1_i; intmat2_i];
    % Discard duplicates
    [Pmat,intmat_i] = uniquerows_integerTable(intmat_full); % Pmat*intmat_i = intmat_full
    intmat_new{i} = intmat_i;
    nterms_arr(i) = size(intmat_i,1);
    % Determine which row in intmat_i (cidcs) relates to which row in the
    % old intmats (ridcs)
    [ridcs,cidcs] = find(Pmat);     % note that Pmat is a permutation matrix
    is_A = ridcs<=nr1;
    ridcs1_arr{i} = ridcs(is_A);    ridcs2_arr{i} = ridcs(~is_A)-nr1;
    cidcs1_arr{i} = cidcs(is_A);    cidcs2_arr{i} = cidcs(~is_A);
end

% Augment the parameter cells with zeros to account for added intervals
Z0 = sparse(size(params1{1},1),size(params1{1},2));
params_new = repmat({Z0},[nterms_arr,1]);   % all cells share the same empty matrix
params1_new = params_new;   params2_new = params_new;
params1_new(cidcs1_arr{:}) = params1(ridcs1_arr{:});
params2_new(cidcs2_arr{:}) = params2(ridcs2_arr{:});

% Declare the output operators
A_out = A_in;
A_out.params = params1_new;
A_out.intmat = intmat_new;
A_out.chkval = true;    % the functional should be properly specified

B_out = B_in;
B_out.params = params2_new;
B_out.intmat = intmat_new;
B_out.chkval = true;

end
