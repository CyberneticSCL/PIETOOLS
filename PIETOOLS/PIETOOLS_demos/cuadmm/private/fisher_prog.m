function [prog,info] = fisher_prog(R,use_bnd,do_solve,slv)
% fisher_prog -- PIESOS_Fisher.m as a FUNCTION, faithful to the shipped script,
% parameterised on the L2-ball radius R and on whether the objective is active.
%
% Fisher: u_t = u_ss + alp*u - bet*u^2, Dirichlet on [0,1], alp=5, bet=-1
% (the shipped values; they reproduce the documented R ~ 4.0479).
% The ball enters ONLY through g = R^2 - <Tx,Tx>, i.e. as a Positivstellensatz
% multiplier -- so bisecting R changes a constant term, not the monomial
% structure.  arXiv 2604.01115 bisects exactly this radius.
%
% use_bnd = false drops gam and piesos_setobj, making the program pure
% FEASIBILITY (c = 0) -- the class cuADMM handles well, and what makes
% bisection on R the right search.
if nargin<2, use_bnd = true;  end
if nargin<3, do_solve = true; end
if nargin<4, slv = 'mosek';   end

clear stateNameGenerator      % name generator carries state between builds
pvar s t
dom = [0,1];  alp = 5;  bet = -1;
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2)+alp*x-bet*x^2;
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];
PIE = convert(PDE);
Top = PIE.T;
f = PIE.f;
x = f.vartab;

d = 1;  pdeg = 4;  eppos = 1;  k = 0;  d_psatz2 = 1;
Zg2 = dmonomials(x,(1:d_psatz2));
g = R^2 - innerprod(Top*x,Top*x);          % the L2 ball

if islogical(use_bnd) && use_bnd
    dpvar gam
    prog = piesos_program(x,gam);
    prog = piesos_setobj(prog,gam);   % OPTIMISATION: c ~= 0
elseif ~islogical(use_bnd)
    gam = use_bnd;                    % FIXED gamma: bound kept, c = 0
    prog = piesos_program(x);
else
    prog = piesos_program(x);         % bound DROPPED entirely (weaker)
end

Zs = cell(d,1);  Pdim_arr = ones(d,1);
for i=1:d, Zs{i} = s.^(0:max(pdeg*i,0)); end
[prog,Pcell] = sosquadvar(prog,Zs,Zs,Pdim_arr,Pdim_arr,'pos');
Pcell{1} = Pcell{1} + eppos*eye(size(Pcell{1}));
Tx = Top*x;
Vx = innerprod(Tx,Tx,Pcell{1});
dV = 2*innerprod(Tx,f,Pcell{1,1});

% psatz term on the derivative condition
lam2_opts.exclude = [1,0,0]';  lam2_opts.deg = pdeg;  lam2_opts.psatz = 0;
[prog,lam2] = piesos_sosvar(prog,Zg2,lam2_opts);
dV_g = dV + lam2*g;

% upper bound on V: kept when the objective is active AND when gamma is
% fixed numerically; dropped only for use_bnd = false.
if ~(islogical(use_bnd) && ~use_bnd)
    V_bnd = gam*innerprod(Tx,Tx) - Vx;
    Z_bnd = polyopvar(f.varname,s,dom);
    Z_bnd.degmat = unique(floor(V_bnd.degmat./2),'rows');
    Q1_opts.deg = pdeg;  Q1_opts.exclude = [1,0,0]';  Q1_opts.psatz = 0:1;
    [prog,W1] = piesos_sosvar(prog,Z_bnd,Q1_opts);
    prog = piesos_eq(prog,V_bnd-W1);
end

% negativity of the derivative
ZQ = polyopvar(f.varname,s,dom);
ZQ.degmat = unique(floor(dV_g.degmat./2),'rows');
Q2_opts.deg = pdeg+1;  Q2_opts.exclude = [1,0,0]';  Q2_opts.psatz = 0:1;
[prog,W2] = piesos_sosvar(prog,ZQ,Q2_opts);
prog = piesos_eq(prog,dV_g+W2);

info = struct('R',R,'use_bnd',use_bnd);
S = sdpshape(prog);
info.m = S.m; info.Kf = S.Kf; info.Ns = S.Ks; info.nblk = S.nblk;
info.c_nnz = nnz(prog.objective);
info.types = strjoin(unique(prog.expr.type),',');

if do_solve
    so.simplify = true;  so.solver = slv;
    t0 = tic;  evalc('prog = lpisolve(prog,so);');  info.t = toc(t0);
    I = prog.solinfo.info;
    Atf=[];bf=[];
    for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
    xv = prog.solinfo.RRx(:);
    info.rel_b = norm(full(Atf'*xv-bf))/norm(full(bf));
    info.normb = norm(full(bf));
    info.pinf=I.pinf; info.dinf=I.dinf; info.numerr=I.numerr; info.feasratio=I.feasratio;
    info.trivial = (abs(info.rel_b-1)<=1e-6) || norm(xv)<=1e-12;
    info.ok = ~I.pinf && ~I.dinf && ~I.numerr && abs(I.feasratio-1)<=0.1 && ~info.trivial;
end
end
