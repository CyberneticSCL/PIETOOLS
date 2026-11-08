function F = horzcat(varargin)
% F = HORZCAT(VARARGIN) takes the horizontal concatenation of distributed
% polynomials listed in VARARGIN, e.g.
%   F = [varargin{1}, varargin{2}, ..., varargin{end}];
%
% INPUTS
% - varargin:   'polyopvar' objects representing distributed polynomials,
%               or numeric arrays, which are converted to 'polyopvar'
%               objects using constant2polyopvar. Empty inputs are
%               ignored. Adjoints (see 'ctranspose') can only be
%               concatenated with other adjoints;
%
% OUTPUTS
% - F:          'polyopvar' object representing the horizontal
%               concatenation of the input functions;
%
% NOTES
% - The output function spaces of the inputs must match. Only scalar or
%   column-vector-valued 'polyopvar' objects are supported, so an error
%   is thrown if the concatenation has more than one column.
% - The horizontal concatenation of adjoints A_i = U_i'*W_i is the adjoint
%   [A_1,...,A_N] = [U_1;...;U_N]'*blkdiag(W_1,...,W_N), so that its
%   product with [V_1;...;V_N] is the sum of the products A_i*V_i. A
%   scalar weight w represents w*eye(n_i), for U_i with n_i rows, and the
%   weight of the concatenation is 1 if all weights are the identity.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - horzcat
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
% DJ, 11/07/2026: Support horzcat for adjoints of polynomials;

% Return an empty 'polyopvar' object if no inputs are specified
if nargin==0
    F = polyopvar();
    return
end
% Adjoints A_i = U_i'*W_i are concatenated as [U_1;...;U_N]'*blkdiag(W_i),
% as [A_1,...,A_N]*[V_1;...;V_N] = sum_i <U_i,W_i*V_i>
if any(cellfun(@(a) isa(a,'polyopvar') && a.is_adjoint,varargin))
    args = varargin(~cellfun(@isempty,varargin));
    if ~all(cellfun(@(a) isa(a,'polyopvar') && a.is_adjoint,args))
        error("polyopvar:adjoint_not_supported","'horzcat' is only supported for adjoint 'polyopvar' objects (u') together with other adjoints, see 'ctranspose'.")
    end
    N = numel(args);
    U = cell(N,1);      W = cell(N,1);      is_id = false(N,1);
    for i=1:N
        if ~args{i}.chkval
            args{i} = initialize(args{i});
        end
        % The polynomial U_i, and the weight W_i (a scalar w as w*eye(n_i))
        ni = size(args{i}.C,1);
        W{i} = args{i}.weight;
        is_id(i) = is_eye(W{i},ni);
        if all(size(W{i})==1) && ni>1
            W{i} = W{i}*eye(ni);
        end
        U{i} = args{i};
        U{i}.is_adjoint = false;        U{i}.weight = 1;
    end
    F = vertcat(U{:});
    F.is_adjoint = true;        F.weight = 1;
    if ~all(is_id)
        F.weight = blkdiag(W{:});
    end
    return
end

% Check that the first input is appropriate
F = varargin{1};
if isempty(F)
    F = horzcat(varargin{2:end});
    return
end
if isnumeric(F)
    F = constant2polyopvar(F);
elseif ~isa(F,'polyopvar')
    error("Concatenation is only supported with objects of type 'polyopvar'.")
end
if ~F.chkval
    F = initialize(F);
end
F.chkval = true;

% Concatenate with the remaining functions
for i=2:nargin
    % Extract the ith function
    B = varargin{i};
    if isempty(B)
        continue
    elseif isnumeric(B)
        B = constant2polyopvar(B);
    elseif ~isa(B,'polyopvar')
        error("Concatenation is only supported with objects of type 'polyopvar'.")
    end
    if ~B.chkval
        B = initialize(B);
    end
    % Express the operators in terms of the same basis
    [F,B] = common_basis(F,B);
    % Make sure the output domains match
    if ~isequal(B.dims.out,F.dims.out)...
        || ~isequal(B.vars.out,F.vars.out)...
            || ~isequal(B.depmat.out,F.depmat.out)
        error("Output function spaces of the polynomials must match for horizontal concatenation.")
    end
    % Concatenate the coefficient operators defining each polynomial
    F.C = [F.C, B.C];
    F.chkval = true;
    if size(F,2)>1  
        error("polyopvar:matrix_valued_not_supported", ...   
            ['Matrix-valued (multi-column) distributed polynomials are ' ...
             'not currently supported -- only scalar or column-vector-' ...
             'valued ''polyopvar'' objects (matdim(2)<=1) are. Horizontal ' ...
             'concatenation of these inputs would produce one with %d ' ...
             'columns.'],size(F,2))                                    
    end                                                                 
end

end