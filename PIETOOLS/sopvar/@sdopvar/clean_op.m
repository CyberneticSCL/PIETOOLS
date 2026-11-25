function [Pop,is_zero] = clean_op(Pop,ztol)
% [POP,IS_ZERO] = CLEAN_OP(POP,ZTOL) takes an 'sdopvar' object POP and
% zeros out coefficients in POP.params below the specified tolerance ZTOL,
% in both the constant part A and the decision variable part B.
%
% INPUTS
% - Pop:    'sdopvar' object representing a PI operator with decision
%           variables;
% - ztol:   (optional) scalar value specifying the tolerance below which
%           coefficients are assumed to be zero. Defaults to 1e-12;
%
% OUTPUTS
% - Pop:    'sdopvar' object representing the same operator as the input
%           up to the specified tolerance, with coefficients in
%           Pop.params.A and Pop.params.B of magnitude at most ztol set
%           to 0;
% - is_zero: true if every coefficient of Pop (across every alpha-index
%           cell of Pop.params) is zero up to ztol, i.e. Pop represents
%           the zero operator for any value of the decision variables;
%           false otherwise.
%
% NOTES
% - The decision variables Pop.Zd are kept, also those that no longer
%   appear in any coefficient.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - clean_op
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 11/24/2026: Initial coding;

% Set default tolerance
if nargin==1
    ztol = 1e-12;
end

% Zero out coefficients below tolerance in A and B of every parameter cell
% (vec C(d) = A + B'*d, sopvar.pdf Sec. 8), and track whether EVERY cell is
% zero up to tolerance (Pop represents the zero operator). Zeroing values
% keeps the canonical multiplier form, so the parameters are written
% directly, without the constructor.
params = Pop.params;
is_zero = true;
for k=1:numel(params.A)
    Ak = params.A{k};   Bk = params.B{k};
    % Test only the stored nonzeros, and rebuild only if some are small:
    % B has a row per decision variable, and abs(B)<=ztol would be true at
    % every zero of B
    if any(abs(nonzeros(Ak))<=ztol)
        Ak = Ak.*(abs(Ak)>ztol);
    end
    if any(abs(nonzeros(Bk))<=ztol)
        Bk = Bk.*(abs(Bk)>ztol);
    end
    params.A{k} = Ak;   params.B{k} = Bk;
    is_zero = is_zero && nnz(Ak)==0 && nnz(Bk)==0;
end
Pop.params = params;

end
