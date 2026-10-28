function prog = lpi_soseq(prog,C,Zd,pos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPI_SOSEQ(PROG,C,ZD,POS) adds to the program PROG the n equations
%
%       C(1,l) + sum_i C(1+i,l) d_i = 0,        l = 1,...,n = size(C,2),
%
% in the decision variables d_i named ZD{i}, as one equality entry
% e = prog.expr.num + 1 of prog.expr: with kept the indices of the nonzero
% columns of C in increasing order, n' their number, N the number of
% decision variables of PROG and pos(i) the row of ZD{i} in
% prog.decvartable,
%
%       type{e} = 'eq',
%       At{e}   = N x n' sparse, At{e}(pos(i),l) = -C(1+i,kept(l)),
%       b{e}    = C(1,kept)', sparse when C is sparse,
%       Z{e}    = sparse(n', numel([prog.vartable; prog.varmat.vartable])),
%
% so that At{e}'*x = b{e}, x the vector of the N decision variables of
% PROG, are the equations of the kept columns.
%
% INPUT
% - prog:  'struct', an SOS or LPI program (sosprogram, lpiprogram) with
%          the fields expr, decvartable, vartable and varmat.vartable;
% - C:     'double', (q+1) x n, sparse or full: row 1 the constant terms,
%          row 1+i the coefficients of ZD{i}, column l equation l;
% - Zd:    'cellstr', the q names of decision variables of PROG;
% - pos:   (optional) 'double', q x 1, pos(i) the first row of
%          prog.decvartable named ZD{i}, used as given; omitted or [], it
%          is computed from the names.
%
% OUTPUT
% - prog:  'struct', PROG with the entry e added to prog.expr.
%
% NOTES
% For a program with symbolic variables (a field symvartable), a C that is
% zero, not a real 'double' or has NaN entries, and, when POS is not given,
% a name repeated in ZD with nonzero rows of C, the entry is that of
% soseq(prog,dpvar(C,zeros(1,0),{},ZD,[1,n])), which is called instead.
% Fewer than three inputs, size(C,1) ~= numel(ZD)+1, numel(POS) ~=
% numel(ZD) and, when POS is not given, a name with a nonzero row of C
% that is not in prog.decvartable are errors.
%
% See also SOSEQ, IMPOSE_EQ_ROWS, COLLECT_EQ_ROWS, LPI_EQ_SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_soseq
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
% Initial coding MMP, 10/02/2026: replaces the soseq call of
%                  'impose_eq_rows'. In the 3-D heatNd build each of its 11
%                  soseq calls spent ~0.54 s in getequation matching the
%                  2.7e6 program names against the expression's (5.9 s of
%                  the 6.86 s getequation total, measured).
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), step comments in the body and the
%                  INPUT and OUTPUT lines of the local function; code
%                  unchanged. The deleted header text was the trace of what
%                  soseq writes for this input, the reason for each case
%                  passed to soseq, the case of names that differ in
%                  trailing blanks, and cost figures against soseq.

% Three inputs at least; ZD a column of q names; C with q+1 rows
if nargin<3
    error("lpi_soseq needs a program, a coefficient matrix and a name list.")
end
Zd = Zd(:);     q = numel(Zd);
n = size(C,2);
if size(C,1)~=q+1
    error("C must have numel(Zd)+1 rows: the constant term, then one row per name.")
end
% The cases where soseq writes something else (see NOTES): soseq itself.
if isfield(prog,'symvartable') || nnz(C)==0 || ~isa(C,'double') || ~isreal(C)
    prog = soseq(prog,dpvar(C,zeros(1,0),{},Zd,[1,n]));
    return
end
% Table rows of the names unless POS is given; a name repeated with nonzero
% rows of C: soseq
if nargin<4 || (isempty(pos) && q>0)
    [pos,dup] = table_rows(prog.decvartable,Zd,C);
    if dup
        prog = soseq(prog,dpvar(C,zeros(1,0),{},Zd,[1,n]));
        return
    end
elseif numel(pos)~=q
    error("pos must give one table row per name in Zd.")
end
pos = pos(:);

% The entry, from the triplets of C. find lists nonzeros only, so the kept
% columns are exactly getequation's any(C,1), and no zero is stored.
[r,c,v] = find(C);
if any(isnan(v))        % soseq's any() ignores NaN: see NOTES
    prog = soseq(prog,dpvar(C,zeros(1,0),{},Zd,[1,n]));
    return
end
r = r(:);   c = c(:);   v = v(:);           % find returns rows for 1 x n C
kept = false(n,1);  kept(c) = true;
newc = cumsum(kept);    nk = nnz(kept);
c = newc(c);                                % kept columns renumbered, in order
isc = (r==1);                               % the constant row
% The new entry e: type, At, b (full when C is full) and Z
e = prog.expr.num + 1;
prog.expr.num = e;
prog.expr.type{e} = 'eq';
prog.expr.At{e} = sparse(pos(r(~isc)-1),c(~isc),-v(~isc),numel(prog.decvartable),nk);
prog.expr.b{e} = sparse(c(isc),1,v(isc),nk,1);
if ~issparse(C)                             % soseq's b is a slice of C
    prog.expr.b{e} = full(prog.expr.b{e});
end
prog.expr.Z{e} = sparse(nk,numel([prog.vartable; prog.varmat.vartable]));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos,dup] = table_rows(dvt,Zd,C)
% POS(i) = first row of the table DVT named Zd{i}; DUP true if a name that C
% uses occurs twice in Zd. Errors on a name C uses that is not in the table,
% as getequation does; names whose rows of C are zero are not needed.
% INPUT  dvt: prog.decvartable, a 'cellstr' or names that 'string'
%             converts;
%        Zd:  q x 1 'cellstr'; C: (q+1) x n 'double'.
% OUTPUT pos: q x 1 'double', the table row of each name, above numel(DVT)
%             for a name not in the table;
%        dup: 'logical'.
% The table names as a column 'cellstr'
if iscellstr(dvt),  known = dvt(:);
else,               known = cellstr(string(dvt(:)));
end
nk = numel(known);
% First occurrences: known precedes Zd, so ia(label) <= nk exactly when the
% name is in the table, and it is then the lowest such row.
[~,ia,ic] = unique([known;Zd]);
pos = ia(ic(nk+1:end));
% The names in use: each in the table, a repeated row flagged
ur = find(any(C,2));    ur = ur(ur>1)-1;    % names with a nonzero row
p = pos(ur);
if any(p>nk)
    error("Decision variable '"+string(Zd{ur(find(p>nk,1))})+"' does not appear in the program.")
end
dup = numel(unique(p))<numel(p);
end
