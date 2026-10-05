function [Pop,is_zero] = clean_op(Pop,ztol)
% [POP,IS_ZERO] = CLEAN_OP(POP,ZTOL) takes a 'sopvar' object POP and zeros
% out coefficients in POP.params below the specified tolerance ZTOL.
%
% INPUTS
% - Pop:    'sopvar' object representing a PI operator;
% - ztol:   scalar value specifying the tolerance below which coefficients
%           are assumed to be zero. Defaults to 1e-12;
%
% OUTPUTS
% - Pop:    'sopvar' object representing the same operator as the input up
%           to the specified tolerance, with coefficients in Pop.params of
%           magnitude smaller than ztol set to 0.
% - is_zero: true if every coefficient of Pop (across every alpha-index
%           cell of Pop.params) is zero up to ztol, i.e. Pop represents the
%           zero operator; false otherwise.
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
% DJ, 10/04/2026: Initial coding;

% Set default tolerance
if nargin==1
    ztol = 1e-12;
end

% Zero out coefficients below tolerance in every parameter cell, and track
% whether EVERY cell is zero up to tolerance (Pop represents the zero
% operator)
params = Pop.params;
is_zero = true;
for k=1:numel(params)
    pk = params{k};
    pk(abs(pk)<=ztol) = 0;
    is_zero = is_zero && nnz(pk)==0;
    params{k} = pk;
end
Pop.params = params;

end
