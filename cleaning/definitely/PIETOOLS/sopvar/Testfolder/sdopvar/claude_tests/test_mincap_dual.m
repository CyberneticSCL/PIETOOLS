function test_mincap_dual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_MINCAP_DUAL runs the pointwise minimum-cap hierarchy (MINCAP_DUAL)
% on the proof program's test cases and on the Volterra slack, and checks
% the two values the proof fixes: the kernel 2 - max(s,t) has cap 1 at lift
% degree 0, and the rank-one target ell_p^* ell_p with deg p = 2 has no
% bounded weight at lift degree 0 (every level unbounded) while its degree-2
% levels stay below the known cap 65700. The other sweeps are printed for
% inspection: the critical Poincare operator G - pi^2 G^2 at D = 1 (no
% bounded weight: the value grows with the input degree and is unbounded
% by degree 12) against D = 2 and 3, and gam - T'T near its optimum.
% Sources: minimum_cap_primal_dual.md Sec. 6, algebraic_hierarchy_duality.md
% Sec. 5, path_b_signed_dual_2026_10_05.md Thm. 3 (Codex outputs, 10/2026).
% Each level is a small SDP (MOSEK; SeDuMi fallback), under a second.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npass = 0;

% (1) K(s,t) = 2 - max(s,t) on [0,1], D = 0 pure: minimum cap 1
P1 = struct('a',0,'b',1,'R0c',[],'R1c',[2;-1],'R2c',[2,-1]);
b1 = sweep('(1) 2 - max(s,t), D = 0 pure (cap 1)',P1,0,0,[2 4 6],[8 16 32]);
if abs(b1(end,end)-1)>2e-2
    error('test_mincap_dual: (1) the cap of 2 - max(s,t) at D = 0 should be 1 (got %.4f).',b1(end,end));
end
fprintf('  passed: (1) cap %.4f at d = 6, M = 32\n',b1(end,end));   npass = npass+1;

% (2) ell_p^* ell_p, p = 30(6t^2 - 6t + 1): D = 0 unbounded at every level; D = 2 below 65700
pc = [30;-180;180];
P2 = struct('a',0,'b',1,'R0c',[],'R1c',pc*pc','R2c',pc*pc');
b2a = sweep('(2) p(s)p(t), D = 0 pure (no bounded weight)',P2,0,0,[2 4],[8 16]);
b2b = sweep('(2) p(s)p(t), D = 2 pure (cap <= 65700)',P2,2,0,[2 4],[16 32]);
if ~all(isinf(b2a(:))) || any(b2b(:)>65700*(1+1e-6))
    error('test_mincap_dual: (2) the rank-one target should be unbounded at D = 0 and below 65700 at D = 2.');
end
fprintf('  passed: (2) D = 0 unbounded at every level; D = 2 levels up to %.4g <= 65700\n',max(b2b(:)));  npass = npass+1;

% (3) critical Poincare G - pi^2 G^2: D = 1 (no bounded weight) against D = 2, 3
pvar s theta
opvar G;    G.I = [0,1];    G.var1 = s;     G.var2 = theta;
G.R.R1 = theta*(1-s);   G.R.R2 = s*(1-theta);
Pc = G - pi^2*(G*G);
P3 = struct('a',0,'b',1,'R0c',pcoef1(Pc.R.R0,s),'R1c',pcoef2(Pc.R.R1,s,theta),'R2c',pcoef2(Pc.R.R2,s,theta));
for D = 1:3
    sweep(sprintf('(3) G - pi^2 G^2, D = %d pure',D),P3,D,0,[4 8 12],[32 128]);
end

% (4) the Volterra slack gam - T'T, T'T kernel 1 - max(s,t), mixed lift
g0 = 4/pi^2;
for del = [1e-1 1e-3]
    P4 = struct('a',0,'b',1,'R0c',g0*(1+del),'R1c',-[1;-1],'R2c',-[1,-1]);
    for D = 0:1
        sweep(sprintf('(4) gam - T''T, gam = 4/pi^2 (1 + %g), D = %d mixed',del,D),P4,D,1,[2 4 8],[16 32]);
    end
end
fprintf('\ntest_mincap_dual passed (%d checks; the sweeps (3) and (4) are for inspection).\n',npass);
end


function B = sweep(lbl,P,D,mult,dl,Ml)
fprintf('\n===== %s\n      M:', lbl);   fprintf(' %10d',Ml);  fprintf('\n');
B = zeros(numel(dl),numel(Ml));
for i = 1:numel(dl)
    fprintf('  d = %2d:',dl(i));
    for j = 1:numel(Ml)
        o = mincap_dual(P,D,dl(i),Ml(j),mult);
        B(i,j) = o.beta;
        flag = ' ';     if ~contains(o.status,'OPTIMAL') && ~isinf(o.beta),  flag = '?';    end
        fprintf(' %9.4g%s',o.beta,flag);
    end
    fprintf('\n');
end
end


function c = pcoef1(p,s)
% coefficient vector of a polynomial in s (index alpha+1 <-> s^alpha)
if isempty(p) || (isa(p,'double') && all(p(:)==0)),  c = [];  return,   end
p = polynomial(p);
if isempty(p.varname), c = double(p); return, end
is = find(strcmp(p.varname,s.varname{1}));
dg = zeros(size(p.degmat,1),1);     if ~isempty(is),    dg = p.degmat(:,is);    end
c = zeros(max(dg)+1,1);
for k = 1:numel(dg),    c(dg(k)+1) = c(dg(k)+1) + p.coefficient(k);    end
end


function C = pcoef2(p,s,th)
% coefficient matrix of a polynomial in (s,theta): (alpha+1,beta+1) <-> s^alpha theta^beta
if isempty(p) || (isa(p,'double') && all(p(:)==0)),  C = [];  return,   end
p = polynomial(p);
if isempty(p.varname), C = double(p); return, end
is = find(strcmp(p.varname,s.varname{1}));  it = find(strcmp(p.varname,th.varname{1}));
nterm = size(p.degmat,1);
da = zeros(nterm,1);    db = zeros(nterm,1);
if ~isempty(is),    da = p.degmat(:,is);    end
if ~isempty(it),    db = p.degmat(:,it);    end
C = zeros(max(da)+1,max(db)+1);
for k = 1:nterm,    C(da(k)+1,db(k)+1) = C(da(k)+1,db(k)+1) + p.coefficient(k);    end
end
