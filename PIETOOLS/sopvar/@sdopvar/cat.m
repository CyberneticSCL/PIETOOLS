function P = cat(dim,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = CAT(DIM,A,B,...) concatenates operators: DIM = 1 is [A; B; ...]
% (vertcat), DIM = 2 is [A, B, ...] (horzcat); any other DIM is an error.
%
% INPUTS
% - dim:        1 or 2;
% - A, B, ...:  operands, as accepted by @sdopvar/vertcat or
%               @sdopvar/horzcat;
%
% OUTPUTS
% - P:          the result of vertcat(A,B,...) or horzcat(A,B,...).
%
% NOTES
% Why: the builtin 'cat' does not call the overloaded vertcat/horzcat, so
% cat(1,A,B), cat(2,A,B) and cat(3,A,B) silently built 2 x 1, 1 x 2 and
% 1 x 1 x 2 ARRAYS of sdopvar objects, which no sdopvar method accepts
% (measured, R2025b). An operator has an output and an input axis only,
% so DIM >= 3 has no operator meaning.
%
% See also VERTCAT, HORZCAT, SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - cat(sdopvar)
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
% Initial coding MMP, 09/29/2026

% Dispatch again, as [A;B] / [A,B] would: the operands decide the class.
if isequal(dim,1)
    P = vertcat(varargin{:});
elseif isequal(dim,2)
    P = horzcat(varargin{:});
else
    error('sdopvar:catDim',['cat(DIM,...) of an sdopvar needs DIM = 1 '...
        '(vertcat) or DIM = 2 (horzcat); any other DIM would build an '...
        'array of operators, which is not an operator.'])
end

end
