function tr = trace_rn_sop(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TR = TRACE_RN_SOP(X) returns the trace of the finite-dimensional part of
% the square PI operator X:
%
%   X 'opvar', 'dopvar'       trace(X.P);
%   X 'opvar2d', 'dopvar2d'   trace(X.R00);
%   X 'copvar', 'cdopvar'     the sum of the traces of the diagonal blocks
%                             X.C{k,k} whose space k has no spatial
%                             variable: for an m x m 'sdopvar' block with
%                             coefficient unvec(A + B'*d), d the decision
%                             variables, sum(A(idx)) + sum(B(:,idx),2)'*d
%                             with idx = (0:m-1)*m + (1:m); for a 'sopvar'
%                             block, trace(params{1}).
%
% INPUT
% - X:    'opvar', 'dopvar', 'opvar2d', 'dopvar2d', 'copvar' or 'cdopvar'
%         object, square.
%
% OUTPUT
% - tr:   for a legacy X, what trace returns for X.P or X.R00 ('double', or
%         'dpvar' for a 'dopvar' or 'dopvar2d' whose part is a 'dpvar');
%         for a container, a 1 x 1 'dpvar' over the decision variables X.Zd
%         (none for a 'copvar'), the constant 0 when every space of X has a
%         spatial variable.
%
% NOTES
% X of another class is an error. For a container, a diagonal block on a
% space without spatial variable that is not square, carries spatial
% variables, or has decision variables other than X.Zd is an error, as is
% an 'sdopvar' block whose A is nonzero and has other than m*m entries.
%
% See also LPI_INEQ, POSLPIVAR_SOP, LPIVAR_CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - trace_rn_sop
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
% Initial coding MMP, 10/06/2026. Library version of the test-folder
%                cx_h2_trace (cx_exec, MMP 09/25/2026), same sums, plus the
%                legacy branches and a check that an R^q diagonal block is
%                square (cx_h2_trace read dims(1) only). cx_h2_trace stays,
%                unchanged, for the frozen 2-D H2 transcriptions.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body; code
%                unchanged. The deleted header text was the use of the
%                output by its callers and a cost note.

% Legacy operators: the trace of the matrix part
if isa(X,'opvar') || isa(X,'dopvar')
    tr = trace(X.P);    return
elseif isa(X,'opvar2d') || isa(X,'dopvar2d')
    tr = trace(X.R00);  return
elseif ~(isa(X,'copvar') || isa(X,'cdopvar'))
    error('trace_rn_sop:class','X must be an opvar, opvar2d, copvar or cdopvar.')
end
% Container: the constant a and the decision coefficients b of the trace,
% over the decision list of X
if isa(X,'cdopvar'),    Zd = X.Zd(:);   else,   Zd = cell(0,1);     end
q = numel(Zd);
a = 0;      b = sparse(q,1);
isR = ~any(X.space_out,2);                  % spaces with no variable: R^q
% Each nonzero diagonal block on a space without spatial variable: check it,
% then add the sum of its diagonal entries (B read on its diagonal columns
% only)
for k = find(isR(:)).'
    blk = X.C{k,k};
    if isempty(blk),    continue,   end
    if ~isempty(blk.vars.out) || ~isempty(blk.vars.in)
        error('trace_rn_sop:vars','Block (%d,%d) of an R^q space carries variables.',k,k)
    end
    m = blk.dims(1);
    if blk.dims(2)~=m
        error('trace_rn_sop:square','Diagonal block (%d,%d) is %d x %d.',k,k,m,blk.dims(2))
    end
    idx = (0:m-1)*m + (1:m);                % diagonal positions in vec
    if isa(blk,'sdopvar')
        A = blk.params.A{1};    B = blk.params.B{1};
        if ~isequal(blk.Zd(:),Zd)
            error('trace_rn_sop:Zd','Block decision list differs from the container''s.')
        end
        if numel(A)==m*m,   a = a + full(sum(A(idx)));
        elseif any(A(:)),   error('trace_rn_sop:A','Nonzero scalar A shorthand in an m>1 block.')
        end                                 % [] or scalar 0: the zero shorthand
        if ~isempty(B), b = b + sum(B(:,idx),2);    end
    else                                    % fixed sopvar block: constants only
        C = blk.params{1};
        if ~isempty(C),  a = a + full(trace(C));   end
    end
end
% The trace a + b'*Zd as a 1 x 1 dpvar
tr = dpvar(sparse([a; b]),zeros(1,0),{},Zd,[1 1]);
end
