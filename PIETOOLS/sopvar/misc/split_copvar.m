function Q = split_copvar(P,rsplit,csplit,opts)
% Q = SPLIT_COPVAR(P,RSPLIT,CSPLIT,OPTS) splits the components of each row
% and column of the 'copvar' object P over several components, returning a
% 'copvar' object Q representing the same operator.
%
% INPUT
% - P:      M x N 'copvar' object (or [], in which case Q=[]);
% - rsplit: (optional) M x 1 cell, with rsplit{i} a cell of index vectors,
%           each specifying the components of the output of row i of P that
%           form a row of Q. Defaults to a single row with all components,
%           for each row i;
% - csplit: (optional) 1 x N cell, with csplit{j} a cell of index vectors,
%           each specifying the components of the input of column j of P
%           that form a column of Q. Defaults to a single column with all
%           components, for each column j;
% - opts:   (optional) 'tensopmat', to return the operators of the blocks
%           as a cell, from which a 'tensopmat' object can be declared;
%
% OUTPUT
% - Q:      'copvar' object with a row for each index vector in rsplit, and
%           a column for each index vector in csplit, e.g. Q.C{r,c} the
%           operator from the components csplit{j}{cc} of column j to the
%           components rsplit{i}{rr} of row i. Where P has a zero block, Q
%           has a zero operator. If opts='tensopmat', Q is instead the
%           cell Q.C of 'sopvar' objects, so that tensopmat(Q) represents
%           the same operator;
%
% NOTES
% - A 'copvar' merges the states on the same variables in a single block,
%   e.g. u_1 and u_2 in v=[D^{a_1}u_1;D^{a_2}u_2]. Splitting v into one
%   component for each state allows them to be treated separately, in the
%   same way as the blocks of a 'tensopmat'. A 'tensopmat' distinguishes
%   each block, even if blocks map to or from the same function space, so
%   that e.g. each state in a distributed polynomial can be a different
%   variable.
% - The components of P in a split component can be in any order.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - split_copvar
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
% DJ, 11/01/2026: Initial coding;

if isempty(P)
    Q = [];
    return
end
if ~isa(P,'copvar')
    error("Operator must be specified as object of type 'copvar'.")
end
as_cell = false;
if nargin>=4
    if ~(ischar(opts) || isstring(opts)) || ~strcmp(opts,'tensopmat')
        error("The option should be specified as 'tensopmat'.")
    end
    as_cell = true;
end
[M,N] = size(P.C);
if nargin<2 || isempty(rsplit)
    rsplit = arrayfun(@(i) {{1:P.dim_out(i)}},(1:M)');
end
if nargin<3 || isempty(csplit)
    csplit = arrayfun(@(j) {{1:P.dim_in(j)}},1:N);
end

% % % Collect the blocks
% Block (rr,cc) of block (i,j) of P is P.C{i,j}(idr,idc), the operator from
% the components idc of the input to the components idr of the output,
% obtained by indexing the sopvar. Row i of P becomes rows r0(i)+(1:nr(i))
% of Q, and column j columns c0(j)+(1:nc(j)). A zero block of P (empty)
% becomes an explicit zero operator on the same spaces, so that every block
% of Q (and of a tensopmat) is declared.
nr = cellfun(@numel,rsplit);        nc = cellfun(@numel,csplit);
ops = cell(sum(nr),sum(nc));
r0 = [0;cumsum(nr(:))];             c0 = [0,cumsum(nc(:))'];
for i=1:M
    for j=1:N
        for rr=1:nr(i)
            for cc=1:nc(j)
                ridcs = rsplit{i}{rr};    cidcs = csplit{j}{cc};
                if ~isempty(P.C{i,j})
                    ops{r0(i)+rr,c0(j)+cc} = P.C{i,j}(ridcs,cidcs);
                    continue
                end
                spc = struct('in',{P.vars(P.space_in(j,:))},'out',{P.vars(P.space_out(i,:))});
                dom = struct('in',P.dom(P.space_in(j,:),:),'out',P.dom(P.space_out(i,:),:));
                ops{r0(i)+rr,c0(j)+cc} = mat2sopvar(zeros(numel(ridcs),numel(cidcs)),spc,dom);
            end
        end
    end
end

% % % Return the blocks as a cell, or declare the 'copvar' object
Q = ops;
if ~as_cell
    Q = copvar(ops);
end

end