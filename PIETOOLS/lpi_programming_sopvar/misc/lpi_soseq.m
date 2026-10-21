function prog = lpi_soseq(prog,C,Zd,pos)
% PROG = LPI_SOSEQ(PROG,C,ZD,POS) imposes the n scalar equations
%       C(1,l) + sum_i C(1+i,l) d_i = 0,        l = 1,...,n = size(C,2),
% on the decision variables d_i named ZD{i}. It appends to prog.expr the
% one 'eq' entry that soseq(prog,dpvar(C,zeros(1,0),{},ZD,[1,n])) appends,
% built from the triplets of C: soseq specialised to the variable-free
% affine equations LPIs impose.
%
% INPUTS
% - prog:   SOS or LPI program structure ('sosprogram', 'lpiprogram');
% - C:      (q+1) x n double, sparse or full; row 1 the constant term, row
%           1+i the coefficients of ZD{i}. Column l is one equation;
% - Zd:     q x 1 cellstr of decision variable names in prog.decvartable;
% - pos:    optional q x 1: pos(i) = the FIRST row of prog.decvartable
%           named Zd{i}, as 'collect_eq_rows' returns it. Trusted - only its
%           length is checked; a wrong pos puts coefficients on the wrong
%           variables. Omitted or [] (q > 0), it is computed here by one
%           unique over the program's names, O((N+q) log(N+q)).
%
% OUTPUT
% - prog:   prog.expr.num + 1 entries; entry e has type{e} = 'eq' and At{e},
%           b{e}, Z{e} equal (isequal, same sparsity) to the ones soseq
%           writes. No other field changes.
%
% NOTES
% For this input soseq -> sosconstr -> getequation (dpvar branch) writes,
% with N = numel(prog.decvartable):
% - the columns of C that are zero (0 = 0) are dropped, the rest kept in
%   order (getequation: col_nums = any(C,1));
% - b = C(1,kept)', a column, sparse when C is;
% - At = N x n' sparse, At(pos(i),l) = -C(1+i,kept(l)), so sossolve's
%   At'*x = b is the equation above;
% - Z = sparse(n', numel([prog.vartable; prog.varmat.vartable])), all
%   zero: one constant monomial per equation.
% These are built here directly. soseq is called instead, as is, when its
% answer is not this: a program with symbolic variables (sosconstr takes
% its symbolic branch); a Zd repeating a name that C uses (combine sums the
% rows, which may cancel a column); an all-zero C (compress returns a full
% coefficient matrix); a C that is not real double (the dpvar constructor
% rejects other classes, getequation conjugates b); and NaN entries (any()
% ignores NaN, so soseq drops a row or column whose only entries are NaN).
% One input class differs: names that differ only in trailing blanks.
% combine pads names into a char array and merges them; here they stay
% distinct. Generated names ('coeff_' + integer) never end in a blank.
%
% Cost, given pos: O(nnz(C) + n) time, about 80 B per nonzero of transient
% memory (triplets, index maps, the output), nothing in N beyond At's
% column pointers. soseq pays a scan of all N names of
% prog.decvartable per call (getequation's ismember), sorts the names of
% the expression (combine) and holds about 9 copies of C.
%
% See also SOSEQ, IMPOSE_EQ_ROWS, COLLECT_EQ_ROWS, LPI_EQ_SDOPVAR.

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
if iscellstr(dvt),  known = dvt(:);
else,               known = cellstr(string(dvt(:)));
end
nk = numel(known);
% First occurrences: known precedes Zd, so ia(label) <= nk exactly when the
% name is in the table, and it is then the lowest such row.
[~,ia,ic] = unique([known;Zd]);
pos = ia(ic(nk+1:end));
ur = find(any(C,2));    ur = ur(ur>1)-1;    % names with a nonzero row
p = pos(ur);
if any(p>nk)
    error("Decision variable '"+string(Zd{ur(find(p>nk,1))})+"' does not appear in the program.")
end
dup = numel(unique(p))<numel(p);
end
