function Psol = getsol_lpivar(prog,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSOL = GETSOL_LPIVAR(PROG,P) returns the decision operator P with
% each of its decision variables replaced by its value in the solved LPI
% program PROG. For an 'sdopvar' operator P with decision variables
% d = P.Zd, whose parameter cell gam stores vec C_gam(d) = A_gam + B_gam'*d
% (A_gam = P.params.A{gam}, B_gam = P.params.B{gam}), the result has the
% parameters
%
%   C_gam = reshape(A_gam + B_gam'*dsol, m, n),
%   m = P.dims(1)*prod(cellfun(@numel,P.ZL)),
%   n = P.dims(2)*prod(cellfun(@numel,P.ZR)),
%
% with dsol(k) the entry of PROG.solinfo.RRx at the first row of
% PROG.decvartable that holds the name P.Zd{k}, and 0 when no row does; a
% name that no row holds keeps its rows of B_gam and remains a decision
% variable. A 'dopvar' container is treated block by block, its 'sopvar'
% and [] blocks kept as they are. A 'sopvar' or 'opvar' P is returned
% unchanged.
%
% INPUT
% - prog:   'struct', a solved LPI program, with the fields
%           decvartable ('cellstr' or string array, the decision variable
%           names) and solinfo.RRx ('double' vector, entry r the value of
%           the r-th name of decvartable);
% - P:      'sdopvar' or 'dopvar' object, the decision operator; or a
%           'sopvar' or 'opvar' object.
%
% OUTPUT
% - Psol:   for an 'sdopvar' P, a 'sopvar' object when PROG.decvartable
%           holds every name of P.Zd, else an 'sdopvar' object whose
%           decision variables are the names it does not hold; for a
%           'dopvar' P, likewise an 'opvar' or a 'dopvar' object, with
%           the spaces, dimensions and domains of P; for the other classes,
%           as stated above.
%
% NOTES
% Errors when PROG has no field decvartable, when solinfo.RRx is empty (no
% solution) or has fewer entries than decvartable, and when P is of
% another class.
%
% See also LPIGETSOL, SUBS_DVAR, SOSGETSOL.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getsol_lpivar
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

% % % Fixed classes: nothing to substitute
if isa(P,'sopvar') || isa(P,'opvar')
    Psol = P;                                   % nothing to substitute
    return
elseif ~isa(P,'sdopvar') && ~isa(P,'dopvar')
    error('getsol_lpivar:badClass',...
        ['No solution extraction for class ''%s''; expected sdopvar, dopvar, '...
         'sopvar or opvar.'],class(P))
end

% % % Solved values, in decvartable order
if ~isstruct(prog) || ~isfield(prog,'decvartable')
    error('getsol_lpivar:badProg','PROG must be an LPI program structure with a decvartable.')
end
if ~isfield(prog,'solinfo') || ~isfield(prog.solinfo,'RRx') || isempty(prog.solinfo.RRx)
    error('getsol_lpivar:notSolved',['No solution seems to have been produced: either '...
        '"lpisolve" has not been called, or the program was found to be infeasible.'])
end
% The names as a cellstr column, and at least one value per name
dtab = prog.decvartable;
if isstring(dtab),  dtab = cellstr(dtab);   end
dtab = dtab(:);
RRx = prog.solinfo.RRx;
if numel(RRx)<numel(dtab)
    error('getsol_lpivar:shortRRx',...
        'solinfo.RRx has %d entries for %d decision variables.',numel(RRx),numel(dtab))
end

% % % One lookup for the whole operator, since every 'sdopvar' block of a
% % % 'dopvar' carries the container list P.Zd: the values d of the names
% % % found, 0 for the others, and the names Zk left as decision variables.
Zd = P.Zd;
if isstring(Zd),    Zd = cellstr(Zd);   end     % names may be strings
Zd = Zd(:);
q = numel(Zd);
% tf(k), loc(k): whether, and at which first row of dtab, the name Zd{k}
% sits: decvartable repeats the name of a Gram entry ((i,j) and (j,i)),
% and the equations of the program use the lowest table row of a name,
% hence the lowest row of each name
if q==0
    tf = false(0,1);    loc = zeros(0,1);
else
    [in,pos] = ismember(dtab,Zd);
    r = find(in);
    % reshape: for a table of one name find returns a ROW, 1x0 when absent,
    % and accumarray rejects row subscripts against a [q,1] size
    loc = accumarray(reshape(pos(r),[],1),reshape(r,[],1),[q,1],@min,0);
    tf = loc>0;
end
d = zeros(q,1);                                 % absent names contribute 0 here
d(tf) = full(RRx(loc(tf)));                     % and keep their rows of B
Zk = Zd(~tf);                                   % names left as decision variables

% An 'sdopvar': substitute in its gamma cells
if isa(P,'sdopvar')
    Psol = subs_block(P,d,tf,Zk);
    return
end

% % % dopvar: one block at a time. A 'sopvar' block and a [] (structurally
% % % zero) block hold no decision variable and are kept as they are.
C = P.C;
for k = 1:numel(C)
    if ~isa(C{k},'sdopvar'),    continue,   end
    if numel(C{k}.Zd)~=q        % O(1) guard; equal lists are a class invariant
        error('getsol_lpivar:badBlockZd',...
            ['Block %d carries %d decision variables, the container %d: the container '...
             'is not on one list (see dopvar/verify).'],k,numel(C{k}.Zd),q)
    end
    C{k} = subs_block(C{k},d,tf,Zk);
end
% The container on the metadata of P: an 'opvar' when every name was found,
% else a 'dopvar' on the names left
meta = metadata(P);                             % validated when P was built
if all(tf)
    Psol = opvar(C,rmfield(meta,'Zd'));         % trusted form: the metadata of P is unchanged
else
    meta.Zd = Zk;                               % every sdopvar block got this array
    Psol = dopvar(C,meta);
end

end


% ========================================================================
% ========================================================================
function Pout = subs_block(P,d,tf,Zk)
% POUT = SUBS_BLOCK(P,D,TF,ZK) returns the 'sdopvar' block P with the
% values D substituted in every gamma cell,
%
%   C_gam = unvec(A_gam + B_gam'*d, m, n),
%
% keeping the rows of B_gam for the names not found (TF false), which are
% ZK.
% INPUT  P:    'sdopvar' object;
%        d:    q x 1 double, the values, 0 for a name not found;
%        tf:   q x 1 logical, true for a name found;
%        Zk:   cellstr, the names not found, in the order of P.Zd.
% OUTPUT Pout: 'sopvar' object when Zk is empty, else an 'sdopvar' object
%              on Zk.
% Size m x n of a coefficient matrix; the parameters of P
q = numel(d);
m = P.dims(1)*prod([cellfun(@numel,P.ZL),1]);  % rows: component outer, ZL inner
n = P.dims(2)*prod([cellfun(@numel,P.ZR),1]);  % columns: component outer, ZR inner
prm = P.params;
A = prm.A;      B = prm.B;
fixed = isempty(Zk);
pC = cell(size(A));                             % fixed parameters, or vec A_gam
pB = cell(size(A));                             % remaining B_gam rows
% Each gamma cell: its constant part a and decision coefficients b
for g = 1:numel(A)
    a = A{g};
    b = [];     if g<=numel(B),     b = B{g};   end
    % [], sparse(0,n) or all zero: no decision content. An all-zero B may be
    % too narrow (q x 1 for a scalar-0 cell), so it is not size-checked;
    % nnz of a sparse B is O(1).
    if isempty(b) || size(b,1)==0 || nnz(b)==0
        v = a;
        bk = sparse(numel(Zk),m*n);
    else
        % Substitute, v = a + b'*d, and keep the rows of b of the names
        % not found
        if size(b,1)~=q || size(b,2)~=m*n
            error('getsol_lpivar:badB',...
                'Cell %d of B is %dx%d; expected %dx%d from Zd, dims, ZL and ZR.',...
                g,size(b,1),size(b,2),q,m*n)
        end
        v = (d.'*b).';                          % (B'd)' = d'B: no transpose of B
        if numel(a)==m*n
            v = v + a;
        elseif ~(isempty(a) || (isscalar(a) && a==0))
            error('getsol_lpivar:badA','Cell %d of A has %d entries; expected %d.',g,numel(a),m*n)
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
        error('getsol_lpivar:badA','Cell %d of A has %d entries; expected %d.',g,numel(v),m*n)
    end
    % The cell: the fixed coefficient matrix, or the new constant part v and
    % the B rows left
    if fixed
        pC{g} = reshape(v,m,n);                 % unvec, column-major
    else
        pC{g} = v;      pB{g} = bk;
    end
end
% The block. Both constructors are handed vars/ZL/ZR/dom unchanged from P,
% which already satisfy the canonical variable order and the canonical
% multiplier form; the substitution only removes content, so neither
% rewrites anything.
if fixed
    Pout = sopvar(pC,P.vars,P.ZL,P.ZR,P.dom,P.dims);
else
    Pout = sdopvar(struct('A',{pC},'B',{pB}),P.vars,Zk,P.ZL,P.ZR,P.dom,P.dims);
end
end
