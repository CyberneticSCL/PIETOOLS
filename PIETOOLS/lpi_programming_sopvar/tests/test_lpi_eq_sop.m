function R = test_lpi_eq_sop()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_LPI_EQ_SOP() checks the dispatch of lpi_eq_sop, and errors at
% the first failed check. On the programs lpiprogram_sop of s on [0,1], of
% (s1,s2) on [0,1] x [0,2] and of (s1,s2,s3) on [0,1] x [0,2] x [0,3]:
%   (1) for 16 cases, lpi_eq_sop(prog,P,opts{:}) returns the program, equal
%       field by field (same_val_sop: 'polynomial' and 'dpvar' objects by
%       their stored arrays), that the routine for the class of P returns
%       when called directly, and adds equality rows; P is a 'cdopvar' from
%       lpivar_cdopvar on L2[s], R x L2[s], L2[s1,s2] or L2[s1,s2,s3], or
%       its sum with its adjoint (lpi_eq_cdopvar); the 'sdopvar' block of
%       the first, of its sum with its adjoint, or of the last
%       (lpi_eq_sdopvar); a 'dopvar' on R x L2[s] or a 'dopvar2d' on
%       L2[s1,s2] from lpivar, or its sum with its adjoint, or the 'dpvar'
%       g - 1 (lpi_eq); opts is {'symmetric'} for the 7 self-adjoint P and
%       {} for the others;
%   (2) for those 7 cases, the call without 'symmetric' adds more rows;
%   (3) the 1x1 'cdopvar' on L2[s] and its 'sdopvar' block add the same
%       rows (At and b);
%   (4) for a fixed 'copvar', a fixed 'sopvar', a 'polynomial', a
%       'double', a 'struct', and the option 'foo' with a 'cdopvar' and
%       with an 'sdopvar', lpi_eq_sop raises the error message and
%       identifier that the routine called directly raises;
%   (5) for a 'dopvar' with the option 'foo', lpi_eq_sop returns the
%       program of lpi_eq without the option.
%
% INPUT
% none.
%
% OUTPUT
% - R:      'struct' with fields n ('double', the number of checks passed)
%           and rows (16x1 'double', the rows added in each case of (1));
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

% No warnings for non-canonical multipliers; variables; check counter
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
pvar s th s1 s2 th1 th2
n = 0;
% % % Programs and decision operators.
p1 = lpiprogram_sop(s,th,[0 1]);
[p1,Pc1] = lpivar_cdopvar(p1,1,{{'s'}},[0 1],2);
[p1,Pd1] = lpivar(p1,[1 1;1 1],[2 2 2]);
[p1,g] = lpidecvar(p1,'g');
[p1,Pc5] = lpivar_cdopvar(p1,[1;1],{{},{'s'}},[0 1],1);   % R^1 x L2[s], 2 x 2
p2 = lpiprogram_sop([s1;s2],[th1;th2],[0 1;0 2]);
[p2,Pc2] = lpivar_cdopvar(p2,1,{{'s1','s2'}},[0 1;0 2],1);
[p2,Pd2] = lpivar(p2,[0 0;0 0;0 0;1 1],1);
v3 = {'s1','s2','s3'};  D3 = [0 1;0 2;0 3];
p3 = lpiprogram_sop(v3,D3);
[p3,Pc3] = lpivar_cdopvar(p3,1,{v3},D3,1);
% A self-adjoint container and its block; the blocks of Pc1 and Pc3
Sc1 = Pc1+Pc1';   Sb1 = Sc1.C{1,1};    Bc1 = Pc1.C{1,1};    Bc3 = Pc3.C{1,1};
% The classes; lpivar on a 2-D program must give a dopvar2d
fprintf('classes: %s %s %s %s %s %s\n',class(Pc1),class(Pd1),class(Pc5),class(Pc2),class(Pd2),class(Pc3));
assert(isa(Pd2,'dopvar2d'),'lpivar on a 2-D program did not give a dopvar2d')

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
  p1, Pd1,          @lpi_eq,         {};
  p1, Pd1+Pd1',     @lpi_eq,         {'symmetric'};
  p2, Pd2,          @lpi_eq,         {};
  p2, Pd2+Pd2',     @lpi_eq,         {'symmetric'};
  p1, g-1,          @lpi_eq,         {}};
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
for k = [2 4 6 8 10 13 15]
    pa = lpi_eq_sop(C{k,1},C{k,2});         % same operator, no option
    assert(nrows(pa)-nrows(C{k,1}) > rows(k),'case %d: symmetric did not reduce the rows (%d vs %d)',...
        k,rows(k),nrows(pa)-nrows(C{k,1}))
    n = n+1;
end
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
  p1, x^2,                                  @lpi_eq,         {};
  p1, 3,                                    @lpi_eq,         {};
  p1, struct('a',1),                        @lpi_eq,         {};
  p1, Pc1,                                  @lpi_eq_cdopvar, {'foo'};
  p1, Pc1.C{1,1},                           @lpi_eq_sdopvar, {'foo'}};
% (4) the same message and identifier as the direct call, which must raise
% an error
for k = 1:size(E,1)
    [p,P,fn,o] = E{k,:};
    ea = err_of(@() lpi_eq_sop(p,P,o{:}));      eb = err_of(@() fn(p,P,o{:}));
    assert(~isempty(eb.message) && strcmp(ea.message,eb.message) && strcmp(ea.identifier,eb.identifier),...
        'error case %d (%s): lpi_eq_sop "%s", direct "%s"',k,class(P),ea.message,eb.message)
    n = n+1;
end
% lpi_eq ignores an unknown option (legacy behaviour), and so does the
% dispatcher: the same program as without it.
assert(same_val_sop(lpi_eq_sop(p1,Pd1,'foo'),lpi_eq(p1,Pd1),'p'),'legacy bad option handled differently')
n = n+1;
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
