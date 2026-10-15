function max_error  = coeff_error(F)
% MAX_ERROR = COEFF_ERROR(F) returns the largest absolute value of the
% coefficients defining the 'polyopvar' object F. If F is the solution of
% the residual f-g of an equality constraint f==g declared with
% `piesos_eq`, this is the largest error in the coefficients.
%
% INPUTS
% - F:      'polyopvar' object, e.g. the output of `piesos_getsol`;
%
% OUTPUTS
% - max_error:  nonnegative scalar, the largest absolute value of the
%               coefficients (kernels of 'intop' operators, or entries of
%               constant coefficients) across all monomials of F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - coeff_error
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
% CRR, 06/09/2026: Initial coding
% DJ, 10/15/2026: Updated to use new polyopvar functions;

% Largest absolute coefficient across the coefficient operators of F
max_error = 0;
for i = 1:numel(F.C.ops)
    max_error = max(max_error,block_coeff_max(F.C.ops{i}));
end

end


%% Local functions
function err = block_coeff_max(op)
% Largest absolute coefficient of the coefficient operator op of a
% 'polyopvar' object: an 'intop' object (sparse coefficients in op.params),
% a constant, or a 'tensopmat' or 'tensopvar' object (nested operators)
err = 0;
if isa(op,'intop')
    for k=1:numel(op.params)
        err = max(err,max_abs(op.params{k}));
    end
elseif isa(op,'tensopmat') || isa(op,'tensopvar')
    for k=1:numel(op.ops)
        err = max(err,block_coeff_max(op.ops{k}));
    end
elseif isa(op,'polynomial') || isa(op,'dpvar')
    err = max_abs(op.C);
elseif isnumeric(op)
    err = max_abs(op);
end
end

function m = max_abs(A)
% Largest absolute value among the nonzero elements of A, or 0
v = nonzeros(abs(A));
m = 0;
if ~isempty(v)
    m = full(max(v));
end
end