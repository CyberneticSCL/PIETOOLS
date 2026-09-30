function Psol = getsol_lpivar_sop(prog,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSOL = GETSOL_LPIVAR_SOP(PROG,P) returns the solved value of a PI operator
% decision variable P of a solved LPI program PROG, for the operator classes
% of both families:
%
%   'sdopvar' -> 'sopvar',   'cdopvar' -> 'copvar'      (new, this file)
%   'dopvar'  -> 'opvar',    'dopvar2d' -> 'opvar2d'    (legacy getsol_lpivar)
%   'sopvar', 'copvar'                                  (returned as is)
%   'opvar', 'opvar2d'                                  (legacy getsol_lpivar)
%
% INPUT
% - prog:   'struct', a solved LPI program (see 'lpiprogram', 'lpisolve').
%           Read: prog.decvartable and prog.solinfo.RRx;
% - P:      PI operator decision variable of one of the classes above.
%
% OUTPUT
% - Psol:   the operator with every decision variable of P that is listed
%           in prog.decvartable replaced by its solved value. When all are
%           listed, which is the case for any variable declared in PROG,
%           Psol is FIXED: 'sopvar' for 'sdopvar', 'copvar' for 'cdopvar'.
%           Names of P absent from prog.decvartable stay decision variables
%           and Psol is then an 'sdopvar'/'cdopvar' on those names only,
%           which is what 'sosgetsol' does for a 'dpvar'.
%
% NOTES
% Values come from prog.solinfo.RRx, whose entries are in decvartable order
% (the order of expr.At rows). prog.solinfo.x is the solver's cone vector:
% it holds a slack per 'ineq' constraint and is in the solver's ordering,
% so it coincides with RRx only for a program with no free variable and no
% inequality.
%
% A decision operator stores each parameter as vec C_gam(d) = A_gam +
% B_gam'*d (sopvar document Sec. 8.1.1), so the solved parameter is
% C_gam = unvec(A_gam + B_gam'*d, m, n), m = dims(1)*prod(|ZL|), n =
% dims(2)*prod(|ZR|), column-major as 'sopvar2sdopvar' vectorizes it. No
% 'dpvar' is formed: legacy 'lpigetsol' on a 'dpvar' of q names goes
% through 'dpvar2poly', measured not to finish in 35 min at q ~ 2.8e5.
%
% COST, q = numel(P.Zd), N = numel(prog.decvartable), nnz = nonzeros of the
% B_gam over all blocks and gamma cells:
% - lookup, ONCE per operator: every 'sdopvar' block of a 'cdopvar' carries
%   the container Zd (class invariant), so one name lookup serves all
%   blocks: decvartable scanned against Zd, O(N log q + q log q) string
%   comparisons, O(N) memory (see dvar_positions);
% - substitution: d'*B_gam per cell, O(nnz(B_gam) + n*m) time. B_gam' is
%   never formed: a transpose of a q-row sparse matrix allocates q+1 column
%   pointers per cell;
% - memory: the q-vector d of values and the lookup's index vectors, plus
%   the result's fixed parameters, which do not depend on q.
% A name that decvartable lists twice (the (i,j) and (j,i) entries of a
% Gram variable) takes its value from its LOWEST row, where getequation
% puts its At row.
%
% See also LPIGETSOL_SOP, SUBS_DVAR_SOP, GETSOL_LPIVAR, LPIGETSOL, SOSGETSOL.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getsol_lpivar_sop
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
% Initial coding MMP, 09/29/2026. Solution extraction for 'sdopvar' and
%                'cdopvar', which had none: every container executive
%                stopped at getsol. The core generalizes 'const_sdopvar' and
%                'recover_dvals' of test_lpi_eq_sdopvar_endtoend.m (sopvar/
%                Testfolder/sdopvar/claude_tests), looking the names up by
%                position instead of through a 'dpvar' of q names. Legacy
%                classes go to 'getsol_lpivar' unchanged, so one executive
%                body can call this for either family.

% % % Legacy and fixed classes: no new code.
if isa(P,'sopvar') || isa(P,'copvar')
    Psol = P;                                   % nothing to substitute
    return
elseif isa(P,'opvar') || isa(P,'dopvar') || isa(P,'opvar2d') || isa(P,'dopvar2d')
    Psol = getsol_lpivar(prog,P);               % legacy, unchanged
    return
elseif ~isa(P,'sdopvar') && ~isa(P,'cdopvar')
    error('getsol_lpivar_sop:badClass',...
        ['No solution extraction for class ''%s''; expected sdopvar, cdopvar, '...
         'sopvar, copvar, opvar, dopvar, opvar2d or dopvar2d.'],class(P))
end

% % % Solved values, in decvartable order. Same test and message as sosgetsol.
if ~isstruct(prog) || ~isfield(prog,'decvartable')
    error('getsol_lpivar_sop:badProg','PROG must be an LPI program structure with a decvartable.')
end
if ~isfield(prog,'solinfo') || ~isfield(prog.solinfo,'RRx') || isempty(prog.solinfo.RRx)
    error('getsol_lpivar_sop:notSolved',['No solution seems to have been produced: either '...
        '"lpisolve" has not been called, or the program was found to be infeasible.'])
end
dtab = prog.decvartable;
if isstring(dtab),  dtab = cellstr(dtab);   end
dtab = dtab(:);
RRx = prog.solinfo.RRx;
if numel(RRx)<numel(dtab)
    error('getsol_lpivar_sop:shortRRx',...
        'solinfo.RRx has %d entries for %d decision variables.',numel(RRx),numel(dtab))
end

% % % One lookup for the whole operator; see COST above.
Zd = P.Zd;
if isstring(Zd),    Zd = cellstr(Zd);   end     % rand_sdopvar names are strings
Zd = Zd(:);
q = numel(Zd);
[tf,loc] = dvar_positions(Zd,dtab);
d = zeros(q,1);                                 % absent names contribute 0 here
d(tf) = full(RRx(loc(tf)));                     % and keep their rows of B
Zk = Zd(~tf);                                   % names left as decision variables

if isa(P,'sdopvar')
    Psol = subs_block(P,d,tf,Zk);
    return
end

% % % cdopvar: one block at a time. A 'sopvar' block and a [] (structurally
% % % zero) block hold no decision variable and are kept as they are.
C = P.C;
for k = 1:numel(C)
    if ~isa(C{k},'sdopvar'),    continue,   end
    if numel(C{k}.Zd)~=q        % O(1) guard; equality of the lists is the class invariant 'verify' checks
        error('getsol_lpivar_sop:badBlockZd',...
            ['Block %d carries %d decision variables, the container %d: the container '...
             'is not on one list (see cdopvar/verify).'],k,numel(C{k}.Zd),q)
    end
    C{k} = subs_block(C{k},d,tf,Zk);
end
meta = metadata(P);                             % validated at P's construction
if all(tf)
    Psol = copvar(C,rmfield(meta,'Zd'));        % trusted form: P's metadata is unchanged
else
    meta.Zd = Zk;                               % every sdopvar block got this array
    Psol = cdopvar(C,meta);
end

end


% ========================================================================
function [tf,loc] = dvar_positions(Zd,dtab)
% loc(k) = LOWEST row of dtab holding Zd{k}, tf(k) = loc(k)>0.
%
% Lowest, because decvartable repeats names: an 'sos' Gram variable of
% size n declares n^2 entries, (i,j) and (j,i) under one name (measured on
% poscopvar: 100 entries, 55 names). getequation puts a name's At row at
% its lowest table row, so RRx there is the value the program's rows
% constrain; sosgetsol (intersect) and ismember(Zd,dtab) read the same.
%
% The table is scanned against Zd, as getequation does, rather than Zd
% looked up in the table: the latter sorts the whole table. Measured,
% N = 3e6 names: 0.65 s against 3.0 s at q = 1e3, 4.06 s against 4.09 s
% at q = 3e6. Time O(N log q + q log q), memory O(N) for tf and loc of
% the scan. Zd must hold distinct names, as every constructor of the
% classes produces; a repeated name would match on its first copy only.
q = numel(Zd);
if q==0
    tf = false(0,1);    loc = zeros(0,1);
    return
end
[in,pos] = ismember(dtab,Zd);
r = find(in);
% reshape: for a table of one name find returns a ROW, 1x0 when absent,
% and accumarray rejects row subscripts against a [q,1] size
loc = accumarray(reshape(pos(r),[],1),reshape(r,[],1),[q,1],@min,0);
tf = loc>0;
end


% ========================================================================
function Pout = subs_block(P,d,tf,Zk)
% Substitute d into every gamma cell of one 'sdopvar' block P,
%
%   C_gam = unvec(A_gam + B_gam'*d, m, n)          (Sec. 8.1.1),
%
% keeping the rows of B_gam for the names not found (tf false), which are
% Zk. Returns a 'sopvar' when Zk is empty, else an 'sdopvar' on Zk. Both
% constructors are handed vars/ZL/ZR/dom unchanged from P, which already
% satisfy the canonical variable order and the canonical multiplier form;
% the substitution only removes content, so neither rewrites anything.
q = numel(d);
m = P.dims(1)*prod([cellfun(@numel,P.ZL),1]);  % rows: component outer, ZL inner
n = P.dims(2)*prod([cellfun(@numel,P.ZR),1]);  % columns: component outer, ZR inner
prm = P.params;
A = prm.A;      B = prm.B;
fixed = isempty(Zk);
pC = cell(size(A));                             % fixed parameters, or vec A_gam
pB = cell(size(A));                             % remaining B_gam rows
for g = 1:numel(A)
    a = A{g};
    b = [];     if g<=numel(B),     b = B{g};   end
    % [], sparse(0,n) or all zero: no decision content. An all-zero B may be
    % too narrow (sopvar2sdopvar of a scalar-0 cell gives q x 1), so it is
    % not size-checked; nnz of a sparse B is O(1).
    if isempty(b) || size(b,1)==0 || nnz(b)==0
        v = a;
        bk = sparse(numel(Zk),m*n);
    else
        if size(b,1)~=q || size(b,2)~=m*n
            error('getsol_lpivar_sop:badB',...
                'Cell %d of B is %dx%d; expected %dx%d from Zd, dims, ZL and ZR.',...
                g,size(b,1),size(b,2),q,m*n)
        end
        v = (d.'*b).';                          % (B'd)' = d'B: no transpose of B
        if numel(a)==m*n
            v = v + a;
        elseif ~(isempty(a) || (isscalar(a) && a==0))
            error('getsol_lpivar_sop:badA','Cell %d of A has %d entries; expected %d.',g,numel(a),m*n)
        end
        bk = [];
        if ~fixed,  bk = b(~tf,:);  end         % O(nnz(b)+q); partial substitution only
    end
    % [] and scalar 0 are the class zero-block shorthand; expand to the
    % stated size so every result parameter has one shape.
    if numel(v)==m*n
        v = sparse(v(:));
    elseif isempty(v) || (isscalar(v) && v==0)
        v = sparse(m*n,1);
    else
        error('getsol_lpivar_sop:badA','Cell %d of A has %d entries; expected %d.',g,numel(v),m*n)
    end
    if fixed
        pC{g} = reshape(v,m,n);                 % unvec, column-major
    else
        pC{g} = v;      pB{g} = bk;
    end
end
if fixed
    Pout = sopvar(pC,P.vars,P.ZL,P.ZR,P.dom,P.dims);
else
    Pout = sdopvar(struct('A',{pC},'B',{pB}),P.vars,Zk,P.ZL,P.ZR,P.dom,P.dims);
end
end
