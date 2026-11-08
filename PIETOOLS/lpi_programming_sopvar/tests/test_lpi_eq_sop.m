function R = test_lpi_eq_sop()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_LPI_EQ_SOP() checks the dispatch of lpi_eq_sop, and errors at
% the first failed check. On the programs lpiprogram of s on [0,1], of
% (s1,s2) on [0,1] x [0,2] and of (s1,s2,s3) on [0,1] x [0,2] x [0,3]:
%   (1) for 12 cases, lpi_eq_sop(prog,P,opts{:}) returns the program, equal
%       field by field (same_val_sop: 'polynomial' and 'dpvar' objects by
%       their stored arrays), that the routine for the class of P returns
%       when called directly, and adds equality rows; P is a 'cdopvar' from
%       lpivar_cdopvar on L2[s], R x L2[s], L2[s1,s2] or L2[s1,s2,s3], or
%       its sum with its adjoint (lpi_eq_cdopvar); the 'sdopvar' block of
%       the first, of its sum with its adjoint, or of the last
%       (lpi_eq_sdopvar); or the 'dpvar' g - 1 (soseq); opts is
%       {'symmetric'} for the 5 self-adjoint P and {} for the others;
%   (2) for those 5 cases, the call without 'symmetric' adds at least as
%       many rows; for the 4 container P that are not self-adjoint (on
%       L2[s] and L2[s1,s2,s3], and their 'sdopvar' blocks), the call with
%       'symmetric' adds fewer rows than without;
%   (3) the 1x1 'cdopvar' on L2[s] and its 'sdopvar' block add the same
%       rows (At and b);
%   (4) for a fixed 'copvar', a fixed 'sopvar', and the option 'foo' with a
%       'cdopvar' and with an 'sdopvar', lpi_eq_sop raises the error
%       message and identifier that the routine called directly raises;
%       for a 'polynomial', a 'double' and a 'struct', the identifier
%       lpi_eq_sop:class;
%   (5) none since 11/08/2026: it checked a 'dopvar' with the option 'foo'
%       against LPI_EQ, both retired.
%
% INPUT
% none.
%
% OUTPUT
% - R:      'struct' with fields n ('double', the number of checks passed)
%           and rows (12x1 'double', the rows added in each case of (1));
%           printed lines with the classes and the counts.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_lpi_eq_sop
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. No header text was deleted. The comment of
%                sopvar2fixed said the A part is added; the code does not
%                add it.
% MMP, 11/01/2026: Check (2) for the container cases: 'collect_eq_rows'
%                now drops equations that are multiples of another, so the
%                call without 'symmetric' on a self-adjoint P drops the
%                adjoint copies as well and may add exactly as many rows.
%                That the option reaches the container routines is checked
%                instead on the 4 container P that are not self-adjoint,
%                where 'symmetric' must add fewer rows. Help: entry (2).
% MMP, 11/04/2026: Programs declared by LPIPROGRAM; LPIPROGRAM_SOP, whose
%                extra input forms (names, containers) only filled fields that
%                container routines do not read, is retired.
% MMP, 11/08/2026: The cases of the opvar classes (09/29/2026) deleted: the
%                'dopvar' and 'dopvar2d' of LPIVAR in (1) and (2), the
%                error cases against LPI_EQ in (4), and (5). The 'dpvar'
%                case is checked against SOSEQ, which lpi_eq_sop calls for
%                it, and a 'polynomial', 'double' or 'struct' P against the
%                identifier lpi_eq_sop:class. LPIVAR, LPI_EQ and the opvar
%                classes are retired (cleaning/private). Help: (1), (2),
%                (4), (5), OUTPUT.

% No warnings for non-canonical multipliers; variables; check counter
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
pvar s th s1 s2 th1 th2
n = 0;
% % % Programs and decision operators.
% p1 = lpiprogram_sop(s,th,[0 1]);                                          % MMP, 11/04/2026 (was)
p1 = lpiprogram(s,th,[0 1]);                                                % MMP, 11/04/2026
[p1,Pc1] = lpivar_cdopvar(p1,1,{{'s'}},[0 1],2);
% [p1,Pd1] = lpivar(p1,[1 1;1 1],[2 2 2]);                                  % MMP, 11/08/2026 (was)
[p1,g] = lpidecvar(p1,'g');
[p1,Pc5] = lpivar_cdopvar(p1,[1;1],{{},{'s'}},[0 1],1);   % R^1 x L2[s], 2 x 2
% p2 = lpiprogram_sop([s1;s2],[th1;th2],[0 1;0 2]);                         % MMP, 11/04/2026 (was)
p2 = lpiprogram([s1;s2],[th1;th2],[0 1;0 2]);                               % MMP, 11/04/2026
[p2,Pc2] = lpivar_cdopvar(p2,1,{{'s1','s2'}},[0 1;0 2],1);
% [p2,Pd2] = lpivar(p2,[0 0;0 0;0 0;1 1],1);                                % MMP, 11/08/2026 (was)
v3 = {'s1','s2','s3'};  D3 = [0 1;0 2;0 3];
% p3 = lpiprogram_sop(v3,D3);                                               % MMP, 11/04/2026 (was)
p3 = lpiprogram(polynomial(v3(:)),[],D3);                                   % MMP, 11/04/2026
[p3,Pc3] = lpivar_cdopvar(p3,1,{v3},D3,1);
% A self-adjoint container and its block; the blocks of Pc1 and Pc3
Sc1 = Pc1+Pc1';   Sb1 = Sc1.C{1,1};    Bc1 = Pc1.C{1,1};    Bc3 = Pc3.C{1,1};
% % The classes; lpivar on a 2-D program must give a dopvar2d               % MMP, 11/08/2026 (was)
% fprintf('classes: %s %s %s %s %s %s\n',class(Pc1),class(Pd1),class(Pc5),class(Pc2),class(Pd2),class(Pc3)); % MMP, 11/08/2026 (was)
% assert(isa(Pd2,'dopvar2d'),'lpivar on a 2-D program did not give a dopvar2d') % MMP, 11/08/2026 (was)
% The classes                                                               % MMP, 11/08/2026
fprintf('classes: %s %s %s %s\n',class(Pc1),class(Pc5),class(Pc2),class(Pc3)); % MMP, 11/08/2026

% % % Dispatch table: {program, P, direct routine, opts or {}}.
C = { ...
  p1, Pc1,          @lpi_eq_cdopvar, {};
  p1, Sc1,          @lpi_eq_cdopvar, {'symmetric'};
  p1, Pc5,          @lpi_eq_cdopvar, {};
  p1, Pc5+Pc5',     @lpi_eq_cdopvar, {'symmetric'};
  p2, Pc2,          @lpi_eq_cdopvar, {};
  p2, Pc2+Pc2',     @lpi_eq_cdopvar, {'symmetric'};
  p3, Pc3,          @lpi_eq_cdopvar, {};
  p3, Pc3+Pc3',     @lpi_eq_cdopvar, {'symmetric'};
  p1, Bc1,          @lpi_eq_sdopvar, {};
  p1, Sb1,          @lpi_eq_sdopvar, {'symmetric'};
  p3, Bc3,          @lpi_eq_sdopvar, {};
  p1, g-1,          @soseq,          {}};                                   % MMP, 11/08/2026
% Rows of the table until 11/08/2026: the opvar cases, retired, and the     % MMP, 11/08/2026
% dpvar case against lpi_eq:                                                % MMP, 11/08/2026
%   p1, Pd1,          @lpi_eq,         {};                                  % MMP, 11/08/2026 (was)
%   p1, Pd1+Pd1',     @lpi_eq,         {'symmetric'};                       % MMP, 11/08/2026 (was)
%   p2, Pd2,          @lpi_eq,         {};                                  % MMP, 11/08/2026 (was)
%   p2, Pd2+Pd2',     @lpi_eq,         {'symmetric'};                       % MMP, 11/08/2026 (was)
%   p1, g-1,          @lpi_eq,         {}};                                 % MMP, 11/08/2026 (was)
% (1) each case against the direct call; the rows it adds
rows = zeros(size(C,1),1);
for k = 1:size(C,1)
    [p,P,fn,o] = C{k,:};
    pa = lpi_eq_sop(p,P,o{:});
    pb = fn(p,P,o{:});
    [tf,where] = same_val_sop(pa,pb,'prog');
    assert(tf,'case %d (%s via %s): differs from the direct call at %s',k,class(P),func2str(fn),where)
    rows(k) = nrows(pa) - nrows(p);
    assert(rows(k)>0,'case %d: no rows added',k)
    n = n+1;
end
% 'symmetric' reaches the routine: the paired cases differ in row count.
% The container writer drops equations that are multiples of another        % MMP, 11/01/2026
% (collect_eq_rows), so without the option the adjoint copies of a          % MMP, 11/01/2026
% self-adjoint P go too and the container counts can tie (cases 2-10);      % MMP, 11/01/2026
% % the legacy writer keeps them (13, 15).                                  % MMP, 11/08/2026 (was)
% for k = [2 4 6 8 10 13 15]                                                % MMP, 11/08/2026 (was)
for k = [2 4 6 8 10]                        % the legacy cases are deleted  % MMP, 11/08/2026
    pa = lpi_eq_sop(C{k,1},C{k,2});         % same operator, no option
%   assert(nrows(pa)-nrows(C{k,1}) > rows(k),'case %d: symmetric did not reduce the rows (%d vs %d)',... % MMP, 11/01/2026 (was)
%       k,rows(k),nrows(pa)-nrows(C{k,1}))                                  % MMP, 11/01/2026 (was)
    more = nrows(pa)-nrows(C{k,1}) > rows(k);                               % MMP, 11/01/2026
    if k<=10,   more = more || nrows(pa)-nrows(C{k,1}) == rows(k);  end     % MMP, 11/01/2026
    assert(more,'case %d: symmetric did not reduce the rows (%d vs %d)',... % MMP, 11/01/2026
        k,rows(k),nrows(pa)-nrows(C{k,1}))                                  % MMP, 11/01/2026
    n = n+1;
end
% The option reaches the container routines: on a P that is not             % MMP, 11/01/2026
% self-adjoint, whose adjoint cells carry other equations, 'symmetric'      % MMP, 11/01/2026
% (one cell of each adjoint pair) adds fewer rows than no option (cases     % MMP, 11/01/2026
% 1, 7, 9, 11: Pc1, Pc3 and their blocks).                                  % MMP, 11/01/2026
for k = [1 7 9 11]                                                          % MMP, 11/01/2026
    pa = lpi_eq_sop(C{k,1},C{k,2},'symmetric');                             % MMP, 11/01/2026
    assert(nrows(pa)-nrows(C{k,1}) < rows(k),'case %d: symmetric did not reach the routine (%d vs %d)',... % MMP, 11/01/2026
        k,nrows(pa)-nrows(C{k,1}),rows(k))                                  % MMP, 11/01/2026
    n = n+1;                                                                % MMP, 11/01/2026
end                                                                         % MMP, 11/01/2026
% A cdopvar and its one block: the same rows through the two routines.
pa = lpi_eq_sop(p1,Pc1);    pb = lpi_eq_sop(p1,Pc1.C{1,1});
[Aa,ba] = rows_of(pa,nrows(p1));    [Ab,bb] = rows_of(pb,nrows(p1));
assert(isequal(Aa,Ab) && isequal(ba,bb),'cdopvar and its block give different rows')
n = n+1;

% % % Errors: each is the routine's own.
pvar x
E = { ...
  p1, copvar({sopvar2fixed(Pc1.C{1,1})}),   @lpi_eq_cdopvar, {};
  p1, sopvar2fixed(Pc1.C{1,1}),             @lpi_eq_sdopvar, {};
  p1, Pc1,                                  @lpi_eq_cdopvar, {'foo'};
  p1, Pc1.C{1,1},                           @lpi_eq_sdopvar, {'foo'}};
% Rows of the table until 11/08/2026, against the retired lpi_eq (checked   % MMP, 11/08/2026
% after the loop by lpi_eq_sop's own identifier):                           % MMP, 11/08/2026
%   p1, x^2,                                  @lpi_eq,         {};          % MMP, 11/08/2026 (was)
%   p1, 3,                                    @lpi_eq,         {};          % MMP, 11/08/2026 (was)
%   p1, struct('a',1),                        @lpi_eq,         {};          % MMP, 11/08/2026 (was)
% (4) the same message and identifier as the direct call, which must raise
% an error
for k = 1:size(E,1)
    [p,P,fn,o] = E{k,:};
    ea = err_of(@() lpi_eq_sop(p,P,o{:}));      eb = err_of(@() fn(p,P,o{:}));
    assert(~isempty(eb.message) && strcmp(ea.message,eb.message) && strcmp(ea.identifier,eb.identifier),...
        'error case %d (%s): lpi_eq_sop "%s", direct "%s"',k,class(P),ea.message,eb.message)
    n = n+1;
end
% A 'polynomial', a 'double' and a 'struct': lpi_eq_sop's own error         % MMP, 11/08/2026
for Q = {x^2, 3, struct('a',1)}                                             % MMP, 11/08/2026
    ea = err_of(@() lpi_eq_sop(p1,Q{1}));                                   % MMP, 11/08/2026
    assert(strcmp(ea.identifier,'lpi_eq_sop:class'),'error case %s: got "%s"',class(Q{1}),ea.identifier) % MMP, 11/08/2026
    n = n+1;                                                                % MMP, 11/08/2026
end                                                                         % MMP, 11/08/2026
% Deleted: (5), lpi_eq_sop against LPI_EQ on a 'dopvar' with the unknown    % MMP, 11/08/2026
% option 'foo' (MMP, 09/29/2026); both are retired.                         % MMP, 11/08/2026
% Summary and output
fprintf('test_lpi_eq_sop: %d checks passed (rows added per case: %s)\n',n,mat2str(rows'));
R = struct('n',n,'rows',rows);
end


% ------------------------------------------------------------------------
function r = nrows(p)
% R = NROWS(P) returns the number of equality rows of the program P, the
% total length of its expr.b entries.
% INPUT  p: 'struct', an LPI program.
% OUTPUT r: 'double'.
r = 0;
for i = 1:p.expr.num,   r = r + numel(p.expr.b{i});     end
end

function [A,b] = rows_of(p,r0)
% [A,B] = ROWS_OF(P,R0) returns the At and b of all expressions of the
% program P, joined, without the first R0 rows (earlier constraints).
% INPUT  p: 'struct', an LPI program; r0: 'double', the number of rows to
%        drop.
% OUTPUT A: 'double' matrix, one column per row kept, the At blocks padded
%        with zero rows to one height; b: 'double' column.
At = [];    bb = [];
% Join the expressions, padding the shorter At with zero rows
for i = 1:p.expr.num
    Ai = p.expr.At{i};
    if ~isempty(At) && size(Ai,1)<size(At,1),   Ai(size(At,1),end) = 0;    end
    if ~isempty(At) && size(At,1)<size(Ai,1),   At(size(Ai,1),end) = 0;    end
    At = [At, Ai];  bb = [bb; p.expr.b{i}];                                 %#ok<AGROW>
end
A = At(:,r0+1:end);     b = bb(r0+1:end);
end

function S = sopvar2fixed(B)
% S = SOPVAR2FIXED(B) returns a fixed, nonzero 'sopvar' on the spaces of
% the 'sdopvar' B: in each parameter cell, the coefficient of the first
% decision variable of B (row 1 of params.B), zero where params.B is
% empty.
% INPUT  B: 'sdopvar' object.
% OUTPUT S: 'sopvar' object.
prm = B.params.A;
for g = 1:numel(prm)
    % Size of the cell: components times monomials, per side
    nL = prod([cellfun(@numel,B.ZL),1]);    nR = prod([cellfun(@numel,B.ZR),1]);
    m = B.dims(1)*nL;   nn = B.dims(2)*nR;
    v = sparse(m*nn,1);
    if ~isempty(B.params.B{g}),     v = v + B.params.B{g}(1,:).';   end
    prm{g} = reshape(v,m,nn);
end
S = sopvar(prm,B.vars,B.ZL,B.ZR,B.dom,B.dims);
end

function e = err_of(f)
% E = ERR_OF(F) returns the message and identifier of the error that F()
% raises (both '' when it raises none).
% INPUT  f: 'function_handle', called with no argument.
% OUTPUT e: 'struct' with fields message and identifier ('char').
e = struct('message','','identifier','');
try,    f();    catch ex,   e = struct('message',ex.message,'identifier',ex.identifier);   end
end
