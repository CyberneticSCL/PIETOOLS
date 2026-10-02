function R = test_lpiprogram_sop(parts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_LPIPROGRAM_SOP(PARTS) tests 'lpiprogram_sop'. PARTS: cellstr
% subset of {'a','b','c','d'} (default all).
%
% (a) Every input form of 'lpiprogram' with <= 2 spatial variables: the
%     output equals lpiprogram's FIELD BY FIELD (polynomial and dpvar
%     objects compared by their stored arrays, not by the elementwise
%     isequal overload); every error case gives lpiprogram's message.
% (b) N = 3, 4: equals the program heatNd_lpi / heatNd_poincare /
%     test_copquadvar_faces built by hand for N > 2 (sosprogram + vartable
%     [vars; vars_dum] + dom), field by field; lpiprogram refuses these.
% (c) The new forms: cellstr names, a copvar / cdopvar / sopvar / sdopvar,
%     and a struct(vars,dom), each equal to the polynomial form with the
%     sorted registry and '_dum' dummies.
% (d) The container workflow on an N-variable program with the legacy
%     lpidecvar, lpi_ineq (scalar), lpisetobj and lpisolve unchanged:
%     - an upper bound on the norm of the N-D Volterra operator
%         (V x)(s) = int_0^{s_1} ... int_0^{s_N} x(t) dt, ||V|| = (2/pi)^N,
%       from  min gam  s.t.  gam*I - V'*V = W, W >= 0 (poscopvar degree 1
%       plus, at N = 1, the product-psatz term as volterra_norm_sop, and
%       at N = 2 the 2N linear face terms), N = 1, 2;
%     - min gam s.t. gam*I - M = W >= 0 for the multiplier M = [2 1;1 2] on
%       L2^2[s1..sN] (gam* = 3), N = 1..4, poscopvar degree 0.
%     Checks: MOSEK numerr 0, rel_b <= 1e-6 on solinfo.RRx rows, Gram
%     blocks PSD; the bound is >= the exact value (soundness), equal to
%     1e-4 where the cone is exact (Volterra N = 1, mult); gam from
%     lpigetsol_sop; the extracted E = gam*I - M - W has coefficient norm
%     between 1 and sqrt(2) times the SDP row residual of its rows.
%
% Initial coding MMP, 09/29/2026
% MMP, 10/01/2026: lpiprogram now takes any number of variables (the cap is
%   removed), so part (a) checks that lpiprogram_sop equals it at N = 3
%   instead of checking that lpiprogram refuses; part (b) still compares
%   N = 3, 4 against the hand-built program, which is now also lpiprogram.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(parts),  parts = {'a','b','c','d'};  end
if ischar(parts),   parts = strsplit(parts,',');    end
R = struct();
if ismember('a',parts),     R.a = part_a();     end
if ismember('b',parts),     R.b = part_b();     end
if ismember('c',parts),     R.c = part_c();     end
if ismember('d',parts),     R.d = part_d();     end
fprintf('test_lpiprogram_sop: ALL PASS (%s)\n',strjoin(parts,','));
end


% ========================================================================
function out = part_a()
pvar s th s1 s2 th1 th2 f1
a = dpvar({'a1','a2'});
PIE = initialize(cx_plant('rd',0.5));
cases = { ...
    {s,th,[0 1]}, ...
    {s,[],[0 1]}, ...
    {s,[0 1]}, ...
    {[s th],[0 1]}, ...                          % n x 2: second column dummies
    {[s1;s2],[th1;th2],[0 1;-1 2]}, ...
    {[s1;s2],[0 1]}, ...                         % one interval for both
    {[s1 th1;s2 th2],[0 1;0 3]}, ...
    {[s1;s2],[],[0 1;0 2]}, ...
    {s,th,[0 1],a}, ...
    {s,th,[0 1],a,f1}, ...
    {s,th,[0 1],f1,a}, ...                       % free and decision swapped
    {s,th,[0 1],f1}, ...                         % a pvar in 4th place: free
    {s,[0 1],a}, ...
    {s,[0 1],a,f1}, ...
    {PIE.vars(:,1),PIE.vars(:,2),PIE.dom}, ...
    {PIE.vars,PIE.dom}};
for k = 1:numel(cases)
    c = cases{k};
    p0 = lpiprogram(c{:});      p1 = lpiprogram_sop(c{:});
    [tf,where] = same_val_sop(p0,p1,'prog');
    assert(tf,'(a) case %d: lpiprogram_sop differs from lpiprogram at %s',k,where)
end
% Non-vacuity: the comparison sees a changed domain and a renamed dummy.
p0 = lpiprogram(s,th,[0 1]);
assert(~same_val_sop(p0,lpiprogram_sop(s,th,[0 2]),'p'),'(a) comparison blind to dom')
assert(~same_val_sop(p0,lpiprogram_sop(s,[],[0 1]),'p'),'(a) comparison blind to vartable')
% Errors: same message as lpiprogram.
bad = { {s}, {s,th,[1 0]}, {[s1;s2],[th1;th2],[0 1;0 1;0 1]}, {s,th,[0 1 2]}, ...
        {s^2,[0 1]}, {s,th,{0,1}}, {s,th,[0 1],3}, {[s1 s2],[0 1;0 1;0 1]}};
for k = 1:numel(bad)
    m0 = err_msg(@() lpiprogram(bad{k}{:}));    m1 = err_msg(@() lpiprogram_sop(bad{k}{:}));
    assert(~isempty(m0) && strcmp(m0,m1),'(a) error case %d: lpiprogram "%s", lpiprogram_sop "%s"',k,m0,m1)
end
% Three variables: lpiprogram refuses, lpiprogram_sop does not.
% (10/01/2026: lpiprogram takes them too; lpiprogram_sop must equal it.)    % MMP, 10/01/2026
pvar s3 th3
% m0 = err_msg(@() lpiprogram([s1;s2;s3],[th1;th2;th3],[0 1]));             % MMP, 10/01/2026 (was)
% assert(contains(m0,'more than 2 spatial variables'),'(a) lpiprogram no longer caps at 2: %s',m0) % MMP, 10/01/2026 (was)
p3 = lpiprogram_sop([s1;s2;s3],[th1;th2;th3],[0 1]);
assert(same_val_sop(lpiprogram([s1;s2;s3],[th1;th2;th3],[0 1]),p3,'p3'), ...
       '(a) 3 variables: lpiprogram_sop differs from lpiprogram')           % MMP, 10/01/2026
assert(isequal(p3.dom,repmat([0 1],3,1)) && isequal(elem_names(p3.vartable),{'s1','s2','s3','th1','th2','th3'}),...
    '(a) 3-variable program wrong')
fprintf('(a) %d forms equal to lpiprogram field by field; %d error messages equal\n',numel(cases),numel(bad));
out = struct('ncases',numel(cases),'nerr',numel(bad));
end


% ========================================================================
function out = part_b()
for N = 3:4
    vars = arrayfun(@(i) sprintf('s%d',i),1:N,'uni',0);
    dom = [zeros(N,1), (1:N)'];
    % The program heatNd_lpi builds for N > 2 (heatNd_lpi.m, 'Program in N
    % variables'), verbatim.
    prog = sosprogram(polynomial([]),dpvar(zeros(0,1)));
    prog.vartable = [prog.vartable; polynomial(vars(:)); polynomial(strcat(vars(:),'_dum'))];
    prog.dom = dom;
    p1 = lpiprogram_sop(polynomial(vars(:)),[],dom);
    [tf,where] = same_val_sop(prog,p1,'prog');
    assert(tf,'(b) N = %d: differs from the hand-built program at %s',N,where)
    p2 = lpiprogram_sop(vars,dom);                  % cellstr registry form
    assert(same_val_sop(prog,p2,'prog'),'(b) N = %d: cellstr form differs',N)
end
fprintf('(b) N = 3, 4 equal to the hand-built N-D program field by field\n');
out = true;
end


% ========================================================================
function out = part_c()
% Registry sorted, domains per variable (distinct, so a row mix-up shows).
vars = {'s1','s2','s3'};    D = [0 1; -1 2; 0 0.5];
ref = lpiprogram_sop(polynomial(vars(:)),[],D);
Ic = eye_copvar_sop(2,{vars},D);                        % copvar
[~,Pc] = lpivar_cdopvar(lpiprogram_sop(vars,D),1,{vars},D,1);  % cdopvar
fm = { {Ic}, {Pc}, {Ic.C{1,1}}, {Pc.C{1,1}}, {struct('vars',{vars},'dom',D)}, {vars,D}, {string(vars),D} };
for k = 1:numel(fm)
    p = lpiprogram_sop(fm{k}{:});
    [tf,where] = same_val_sop(ref,p,'prog');
    assert(tf,'(c) form %d (%s) differs at %s',k,class(fm{k}{1}),where)
end
% A container over two spaces with different variables: registry order.
Mc = mat2copvar_sop(ones(2,2),struct('out',[1;1],'in',[1;1]),...
    struct('out',{{{'b'},{'a'}}},'in',{{{'b'},{'a'}}}),[0 3;0 1]);
p = lpiprogram_sop(Mc);
assert(isequal(elem_names(p.vartable),{'a','b','a_dum','b_dum'}) && isequal(p.dom,[0 3;0 1]),...
    '(c) registry order or domains wrong: %s, %s',strjoin(elem_names(p.vartable),','),mat2str(p.dom))
% With decision variables after a container.
p = lpiprogram_sop(Ic,dpvar({'g'}));
assert(isequal(p.decvartable,{'g'}),'(c) decision variables after a container lost')
% Non-vacuity: a different domain is seen.
assert(~same_val_sop(ref,lpiprogram_sop(vars,[0 1;0 2;0 0.5]),'p'),'(c) comparison blind to dom')
fprintf('(c) %d new input forms equal to the polynomial form\n',numel(fm));
out = true;
end


% ========================================================================
function out = part_d()
% Two N-variable LPIs, solved with MOSEK through the legacy lpidecvar,
% lpi_ineq (scalar), lpisetobj and lpisolve on an lpiprogram_sop program:
%  'volterra' (N = 1, 2): min gam s.t. gam*I - V'V = W, W >= 0 from
%     poscopvar degree 1 plus psatz terms (below); ||V|| = (2/pi)^N.
%     Degree 1 at N = 3 is 2 Grams of 1000 (2.0e6 variables, m 72131,
%     build 105 s, measured), so N <= 2 here;
%  'mult' (N = 1..4): min gam s.t. gam*I - M = W, W >= 0 from poscopvar
%     degree 0, M the constant multiplier [2 1; 1 2] on L2^2[s1..sN]:
%     gam* = 3, the largest eigenvalue of [2 1; 1 2].
sopts = struct('solver','mosek','simplify',false);
out = struct('prob',{},'N',{},'gam',{},'exact',{},'rel_b',{},'ce_rr',{},'ndv',{},'m',{},'Ks',{},'t',{});
C = [repmat({'volterra'},2,1), {1;2}; repmat({'mult'},4,1), {1;2;3;4}];
for c = 1:size(C,1)
    [pb,N] = C{c,:};
    t0 = tic;
    vars = arrayfun(@(i) sprintf('s%d',i),1:N,'uni',0);
    D = repmat([0 1],N,1);
    if strcmp(pb,'volterra')
        % V: one lower-integral cell (gamma = (2,...,2)), constant kernel 1.
        prm = repmat({sparse(1,1)},[3*ones(1,N),1]);
        g2 = num2cell(2*ones(1,N));     prm{sub2ind([3*ones(1,N),1],g2{:})} = sparse(1);
        V = copvar({sopvar(prm,struct('in',{vars},'out',{vars}),repmat({0},1,N),repmat({0},1,N),...
            struct('in',D,'out',D),[1 1])});
        M = V'*V;   nc = 1;     deg = 1;    ex = (2/pi)^N;
    else
        M = mat2copvar_sop([2 1;1 2],2,{vars},D);
        nc = 2;     deg = 0;    ex = 3;
    end
    prog = lpiprogram_sop(M);                           % container form
    [prog,gam] = lpidecvar(prog,'gam');
    prog = lpi_ineq(prog,gam);                          % legacy scalar path
    [prog,W] = poscopvar(prog,nc,{vars},D,deg);
    if strcmp(pb,'volterra') && N==1
        [prog,W1] = poscopvar(prog,nc,{vars},D,deg,struct('psatz',1));
        W = W + W1;
    elseif strcmp(pb,'volterra')
        % The product generator certifies nothing here (MOSEK numerr 2,
        % rel_b 1.2e-4, measured N = 2); the 2N linear face generators
        % (copquadvar psatz = 2i+1 / 2i+2, as heatNd_lpi) instead.
        for i = 1:N
            for e = 1:2
                [prog,W1] = poscopvar(prog,nc,{vars},D,deg,struct('psatz',2*i+e));
                W = W + W1;
            end
        end
    end
    E = gam - M - W;                                    % dpvar - copvar: gam*I - M
    ne0 = prog.expr.num;
    prog = lpi_eq_sop(prog,E,'symmetric');
    rows = ne0+1:prog.expr.num;
    prog = lpisetobj(prog,gam);
    S0 = cx_shape(prog);
    assert(S0.ndv<=1e5,'(d) %s N=%d: %d decision variables, above this test''s budget',pb,N,S0.ndv)
    evalc('prog = lpisolve(prog,sopts);');
    info = prog.solinfo.info;
    [rb,~,prel] = cx_resid(prog);
    g = double(lpigetsol_sop(prog,gam));
    Es = lpigetsol_sop(prog,E);     Ws = lpigetsol_sop(prog,W);
    % Every row of the constraint is a canonical coefficient of E, one per
    % adjoint pair ('symmetric'), so ||r_rows|| <= ||coef E(sol)|| <=
    % sqrt(2) ||r_rows||: extraction reproduces the solver's row residual.
    rr = row_res(prog,rows);    ce = coef_norm(Es);
    S = cx_shape(prog);
    if strcmp(pb,'volterra'),   val = sqrt(g);  else,   val = g;    end
    out(c) = struct('prob',pb,'N',N,'gam',val,'exact',ex,'rel_b',rb,'ce_rr',ce/max(rr,realmin),...
        'ndv',S.ndv,'m',S.m,'Ks',S.Ks,'t',toc(t0));
    fprintf(['(d) %-8s N=%d: bound %.8f, exact %.8f (ratio %.6f); MOSEK numerr %d pinf %d, rel_b %.2e, '...
             'psd_relmin %.1e; ||coef E(sol)|| %.2e / ||rows|| %.2e = %.4f; ndv %d, m %d, Ks %s; %.1f s\n'],...
        pb,N,val,ex,val/ex,info.numerr,info.pinf,rb,prel,ce,rr,ce/max(rr,realmin),S.ndv,S.m,mat2str(S.Ks),toc(t0));
    assert(info.numerr==0 && info.pinf==0 && rb<=1e-6 && prel>=-1e-8,'(d) %s N=%d: solve not certified',pb,N)
    assert(val>=ex*(1-1e-6),'(d) %s N=%d: bound %.8g below the exact value %.8g',pb,N,val,ex)
    if strcmp(pb,'mult') || N==1
        assert(val<=ex*(1+1e-4),'(d) %s N=%d: bound %.8g not within 1e-4 of %.8g',pb,N,val,ex)
    end
    assert(isa(Es,'copvar') && isa(Ws,'copvar'),'(d) %s N=%d: extracted operators not fixed',pb,N)
    assert(ce>=rr*(1-1e-6) && ce<=sqrt(2)*rr*(1+1e-6),...
        '(d) %s N=%d: ||coef E(sol)|| = %.3g does not match the SDP row residual %.3g',pb,N,ce,rr)
end
end


% ========================================================================
function nm = elem_names(p)
% Variable name of each element of a column of pvars, in element order
% (p.varname is the sorted list of the array's variables).
nm = cell(1,numel(p));
for i = 1:numel(p),     v = p(i).varname;   nm{i} = v{1};     end
end


% ========================================================================
function r = row_res(prog,rows)
% ||At_i' x - b_i|| over the expressions ROWS, x = solinfo.RRx (decvartable
% order, the order of At rows).
x = prog.solinfo.RRx(:);    r2 = 0;
for i = rows
    At = prog.expr.At{i};   b = prog.expr.b{i};
    r2 = r2 + sum((At'*x(1:size(At,1)) - b).^2);
end
r = sqrt(full(r2));
end

function c = coef_norm(P)
c = 0;
for k = 1:numel(P.C)
    if isempty(P.C{k}),     continue,   end
    prm = P.C{k}.params;
    for g = 1:numel(prm),   c = c + sum(abs(prm{g}(:)).^2);    end
end
c = sqrt(full(c));
end

function m = err_msg(f)
m = '';
try,    f();    catch e,    m = e.message;  end
end
