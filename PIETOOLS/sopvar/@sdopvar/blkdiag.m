function Pblkdiag = blkdiag(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pblkdiag = blkdiag(A,B,...) places sdopvar operators on the block diagonal,
% so that Pblkdiag maps the stacked input (x_A;x_B;...) to the stacked
% output (A*x_A; B*x_B; ...).
%
% INPUT
% varargin: sdopvar objects sharing vars and dom. Dimensions need not match.
%           A fixed 'sopvar' operand is accepted and promoted;
%
% OUTPUT
% Pblkdiag: sdopvar object of dimension [sum of output dims, sum of input dims]
%
% NOTES:
% Built as the block matrix whose diagonal holds the operands and whose
% off-diagonal blocks are zero operators, using @sdopvar/horzcat and
% @sdopvar/vertcat. Those already reconcile the monomial bases and decision
% variable lists of all N operands in one pass, so nothing is gained by
% doing the index arithmetic again here, and the off-diagonal blocks cost
% almost nothing: an all-zero parameter is a sparse allocation with no
% stored entries.
%
% MATLAB's builtin blkdiag cannot be used on these objects: it interleaves
% numeric zero matrices among the operands and hands them to horzcat, which
% has no meaning for an operator on L_2.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/07/2026

% Deal with single input case
if nargin==1
    Pblkdiag = varargin{1};
    return
end

% A fixed 'sopvar' operand is promoted to a decision operator with a zero B,
% taking its decision variable list from whichever operand already has one.
iss = cellfun(@(x) isa(x,'sopvar'),varargin);
if any(iss)
    isd = find(cellfun(@(x) isa(x,'sdopvar'),varargin),1);
    if isempty(isd),    Zd_p = cell(0,1);
    else,               Zd_p = varargin{isd}.Zd;
    end
    for ii = find(iss)
        varargin{ii} = sopvar2sdopvar(varargin{ii},Zd_p);
    end
end
ops = varargin;
N = numel(ops);

if ~all(cellfun(@(x) isa(x,'sdopvar'),ops))
    error('blkdiag is currently supported only for sdopvar and sopvar objects');
end

% Row i holds zeros except for operand i on the diagonal. horzcat and
% vertcat perform the compatibility checks on vars and dom.
mk = cellfun(@(P) P.dims(1),ops);
nk = cellfun(@(P) P.dims(2),ops);
rows = cell(1,N);
blocks = cell(1,N);
for i = 1:N
    for j = 1:N
        if i==j
            blocks{j} = ops{i};
        else
            blocks{j} = zero_like(ops{i},[mk(i),nk(j)]);
        end
    end
    rows{i} = horzcat(blocks{:});
end
Pblkdiag = vertcat(rows{:});

end


%%
function Z = zero_like(P,dims)
% Zero operator of the requested dimension, on the variables, domain,
% monomial bases and decision variable list of P. The bases are inherited so
% that the concatenations see one common basis wherever possible; an all-zero
% coefficient is correct on any basis.

NL = prod([cellfun(@numel,P.ZL),1]);
NR = prod([cellfun(@numel,P.ZR),1]);
nC = dims(1)*NL*dims(2)*NR;
q  = numel(P.Zd);

params.A = cell(size(P.params.A));
params.B = cell(size(P.params.B));
for k = 1:numel(params.A)
    % 'sparse' allocates only the column pointers, so a zero block costs
    % nothing even when there are millions of decision variables.
    params.A{k} = sparse(nC,1);
    params.B{k} = sparse(q,nC);
end

Z = sdopvar(params,P.vars,P.Zd,P.ZL,P.ZR,P.dom,dims);

end
