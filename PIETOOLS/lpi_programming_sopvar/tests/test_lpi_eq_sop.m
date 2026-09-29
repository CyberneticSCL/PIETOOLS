function R = test_lpi_eq_sop()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_LPI_EQ_SOP() tests the dispatch of 'lpi_eq_sop': for every class
% it routes, the program it returns equals, FIELD BY FIELD, the program the
% routine of its dispatch table returns when called directly, with and
% without the 'symmetric' option; errors are the routine's own.
%
% Classes: cdopvar (1-D; 2-D; 3 variables on an lpiprogram_sop program; an
% R^1 x L2 2 x 2 container), sdopvar (a block), dopvar (lpivar, 1-D),
% dopvar2d (lpivar, 2-D), dpvar (soseq). Non-vacuity: every call adds
% rows; 'symmetric' changes the row count, so the option is passed on; the
% cdopvar and its lone block give the same rows through the two routines.
% Refused: copvar, sopvar (fixed), polynomial, double, struct; a bad option
% for each family.
%
% Initial coding MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
Sc1 = Pc1+Pc1';   Sb1 = Sc1.C{1,1};    Bc1 = Pc1.C{1,1};    Bc3 = Pc3.C{1,1};
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
fprintf('test_lpi_eq_sop: %d checks passed (rows added per case: %s)\n',n,mat2str(rows'));
R = struct('n',n,'rows',rows);
end


% ------------------------------------------------------------------------
function r = nrows(p)
r = 0;
for i = 1:p.expr.num,   r = r + numel(p.expr.b{i});     end
end

function [A,b] = rows_of(p,r0)
% At and b of all expressions, joined; the first r0 rows (earlier
% constraints) dropped.
At = [];    bb = [];
for i = 1:p.expr.num
    Ai = p.expr.At{i};
    if ~isempty(At) && size(Ai,1)<size(At,1),   Ai(size(At,1),end) = 0;    end
    if ~isempty(At) && size(At,1)<size(Ai,1),   At(size(Ai,1),end) = 0;    end
    At = [At, Ai];  bb = [bb; p.expr.b{i}];                                 %#ok<AGROW>
end
A = At(:,r0+1:end);     b = bb(r0+1:end);
end

function S = sopvar2fixed(B)
% The block at d = 0: a fixed, nonzero 'sopvar' (A part plus a row of B).
prm = B.params.A;
for g = 1:numel(prm)
    nL = prod([cellfun(@numel,B.ZL),1]);    nR = prod([cellfun(@numel,B.ZR),1]);
    m = B.dims(1)*nL;   nn = B.dims(2)*nR;
    v = sparse(m*nn,1);
    if ~isempty(B.params.B{g}),     v = v + B.params.B{g}(1,:).';   end
    prm{g} = reshape(v,m,nn);
end
S = sopvar(prm,B.vars,B.ZL,B.ZR,B.dom,B.dims);
end

function e = err_of(f)
e = struct('message','','identifier','');
try,    f();    catch ex,   e = struct('message',ex.message,'identifier',ex.identifier);   end
end
