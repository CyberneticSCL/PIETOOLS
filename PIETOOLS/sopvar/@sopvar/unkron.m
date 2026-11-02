function Xs = unkron(X,n)
% Xs = UNKRON(X,N) extracts the 'sopvar' object Xs such that X is the
% Kronecker product kron(Xs,speye(N)), if X has this form, and throws an
% error otherwise. In kron(Xs,speye(N)) the same operator Xs acts on N
% separate copies of the components, with the components of the copies
% (index i=1,...,N) interleaved, so that the output component (a-1)*N+i is
% the component a of Xs applied to the components (b-1)*N+i, b=1,...,s of
% the input, see also kron.
%
% INPUTS
% - X:      r*N x s*N 'sopvar' object, from L_2^{s*N}[S1,S3] to
%           L_2^{r*N}[S2,S3];
% - N:      positive integer, the number of copies of the components;
%
% OUTPUTS
% - Xs:     r x s 'sopvar' object, from L_2^s[S1,S3] to L_2^r[S2,S3], such
%           that X=kron(Xs,speye(N)), equal to the block of X with the
%           components (a-1)*N+1 of the output and (b-1)*N+1 of the input;
%
% NOTES
% - Xs is read from the components with i=1, and X is then required to
%   coincide with kron(Xs,speye(N)), up to a tolerance of 1e-10 relative to
%   the largest coefficient of each parameter of X. Otherwise, the error
%   'sopvar:unkron:notKron' is thrown.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - unkron (sopvar)
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
% DJ, 11/01/2026: Initial coding.

% % % Check the inputs
if nargin~=2
    error("Exactly two inputs are required.")
end
if ~isa(X,'sopvar')
    error("Operator must be specified as object of type 'sopvar'.")
end
if ~isnumeric(n) || ~isscalar(n) || n<1 || n~=round(n)
    error("The number of copies must be specified as a positive integer.")
end
if any(mod(X.dims,n)~=0)
    error('sopvar:unkron:notKron',"The dimensions of the operator are not multiples of %d.",n)
end

% % % Extract the block with the first copy
% The rows of the parameters are (component, monomial), monomial fastest,
% so the rows of the components (a-1)*n+1 are those in the selected
% blocks of nnZL_arr consecutive rows, and likewise for the columns
nnZL_arr = prod(cellfun(@numel,X.ZL));
nnZR_arr = prod(cellfun(@numel,X.ZR));
dims = X.dims/n;
ridcs = reshape((0:n:X.dims(1)-1)*nnZL_arr+(1:nnZL_arr)',[],1);
cidcs = reshape((0:n:X.dims(2)-1)*nnZR_arr+(1:nnZR_arr)',[],1);
params_new = X.params;
for lidx=1:numel(params_new)
    if isempty(params_new{lidx}) || (isscalar(params_new{lidx}) && params_new{lidx}==0)
        continue                            % zero block stays zero
    end
    params_new{lidx} = params_new{lidx}(ridcs,cidcs);
end
Xs = sopvar(params_new,X.vars,X.ZL,X.ZR,X.dom,dims);

% % % Check that X is the Kronecker product of Xs with the identity
X_kron = kron(Xs,speye(n));
param_size = [X.dims(1)*nnZL_arr,X.dims(2)*nnZR_arr];
for lidx=1:numel(X.params)
    param_X = X.params{lidx};           param_kron = X_kron.params{lidx};
    % Expand a zero block in the shorthand [] or 0
    if ~isequal(size(param_X),param_size)
        param_X = sparse(param_size(1),param_size(2));
    end
    if ~isequal(size(param_kron),param_size)
        param_kron = sparse(param_size(1),param_size(2));
    end
    err_max = max([0;abs(nonzeros(param_X-param_kron))]);
    if err_max>1e-10*max(1,max([0;abs(nonzeros(param_X))]))
        error('sopvar:unkron:notKron',"The operator is not of the form kron(Xs,speye(%d)).",n)
    end
end

end
