function Psol = getsol_lpivar_sop(prog,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSOL = GETSOL_LPIVAR_SOP(PROG,P) returns the decision operator P with
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
% variable. A 'cdopvar' container is treated block by block, its 'sopvar'
% and [] blocks kept as they are. A 'sopvar' or 'copvar' P is returned
% unchanged; an 'opvar', 'dopvar', 'opvar2d' or 'dopvar2d' P is passed to
% GETSOL_LPIVAR.
%
% INPUT
% - prog:   'struct', an LPI program solved by LPISOLVE, with the fields
%           decvartable ('cellstr' or string array, the decision variable
%           names) and solinfo.RRx ('double' vector, entry r the value of
%           the r-th name of decvartable);
% - P:      'sdopvar' or 'cdopvar' object, the decision operator; or a
%           'sopvar', 'copvar', 'opvar', 'dopvar', 'opvar2d' or 'dopvar2d'
%           object.
%
% OUTPUT
% - Psol:   for an 'sdopvar' P, a 'sopvar' object when PROG.decvartable
%           holds every name of P.Zd, else an 'sdopvar' object whose
%           decision variables are the names it does not hold; for a
%           'cdopvar' P, likewise a 'copvar' or a 'cdopvar' object, with
%           the spaces, dimensions and domains of P; for the other classes,
%           as stated above.
%
% NOTES
% Errors when PROG has no field decvartable, when solinfo.RRx is empty (no
% solution) or has fewer entries than decvartable, and when P is of
% another class.
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and headers
%                of the local functions; code unchanged. The deleted
%                header text was the class table, the reason solinfo.RRx
%                is read rather than solinfo.x, the cost analysis with a
%                measured run time, and the design rationale and
%                measurements of the local-function headers.

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
% The names as a cellstr column, and at least one value per name
dtab = prog.decvartable;
if isstring(dtab),  dtab = cellstr(dtab);   end
dtab = dtab(:);
RRx = prog.solinfo.RRx;
if numel(RRx)<numel(dtab)
    error('getsol_lpivar_sop:shortRRx',...
        'solinfo.RRx has %d entries for %d decision variables.',numel(RRx),numel(dtab))
end

% % % One lookup for the whole operator, since every 'sdopvar' block of a
% % % 'cdopvar' carries the container list P.Zd: the values d of the names
% % % found, 0 for the others, and the names Zk left as decision variables.
Zd = P.Zd;
if isstring(Zd),    Zd = cellstr(Zd);   end     % rand_sdopvar names are strings
Zd = Zd(:);
q = numel(Zd);
[tf,loc] = dvar_positions(Zd,dtab);
d = zeros(q,1);                                 % absent names contribute 0 here
d(tf) = full(RRx(loc(tf)));                     % and keep their rows of B
Zk = Zd(~tf);                                   % names left as decision variables

% An 'sdopvar': substitute in its gamma cells
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
% The container on the metadata of P: a 'copvar' when every name was found,
% else a 'cdopvar' on the names left
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
% [TF,LOC] = DVAR_POSITIONS(ZD,DTAB) returns, for each name Zd{k}, the
% first row LOC(k) of DTAB that holds it (0 for none) and TF(k) = LOC(k)>0.
% INPUT  Zd:   q x 1 cellstr of distinct names (a repeated name is found on
%              its first copy only);
%        dtab: N x 1 cellstr, the decision variable table.
% OUTPUT tf:   q x 1 logical;
%        loc:  q x 1 double.
% No names: empty outputs
q = numel(Zd);
if q==0
    tf = false(0,1);    loc = zeros(0,1);
    return
end
% Each table row against Zd, then the lowest row of each name: decvartable
% repeats the name of a Gram entry ((i,j) and (j,i)), and getequation puts
% the At row of a name at its lowest table row
[in,pos] = ismember(dtab,Zd);
r = find(in);
% reshape: for a table of one name find returns a ROW, 1x0 when absent,
% and accumarray rejects row subscripts against a [q,1] size
loc = accumarray(reshape(pos(r),[],1),reshape(r,[],1),[q,1],@min,0);
tf = loc>0;
end


% ========================================================================
function Pout = subs_block(P,d,tf,Zk)
% POUT = SUBS_BLOCK(P,D,TF,ZK) returns the 'sdopvar' block P with the
% values D substituted in every gamma cell,
%
%   C_gam = unvec(A_gam + B_gam'*d, m, n)          (Sec. 8.1.1),
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
    % too narrow (sopvar2sdopvar of a scalar-0 cell gives q x 1), so it is
    % not size-checked; nnz of a sparse B is O(1).
    if isempty(b) || size(b,1)==0 || nnz(b)==0
        v = a;
        bk = sparse(numel(Zk),m*n);
    else
        % Substitute, v = a + b'*d, and keep the rows of b of the names
        % not found
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
