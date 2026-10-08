function test_lpi_ineq_sop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_LPI_INEQ_SOP checks the container inequality 'lpi_ineq_sop'.
% (1) A fixed positive operator, V'*V + L'*L on [0,1] (lower kernel
%     2 - s), is certified with the default sizing: the solved Pop
%     reproduces it (rel_b at solver accuracy), the reader chose lift 2.
% (2) The 1-D heat stability LPI of Cor. 35 (heatNd_pie, P the tensor
%     basis of degree 1) posed with two inequalities: feasible at kappa = 9
%     and infeasible at kappa = 11 (kappa* = pi^2), the reader's degrees
%     being R (2,1) and Q (2,2) with the plain and product terms.
% (3) 2-D: T'*T of the heat PIE is certified with the default face terms.
% (4) Options: 'none', 'faces', a given 'deg', prune off: all build, with
%     the expected decision counts ordered.
% (5) A legacy 'dopvar' goes to 'lpi_ineq' unchanged.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
sopts = struct('solver','mosek','simplify',false);
if isempty(which('mosekopt')),  sopts.solver = 'sedumi';   end
npass = 0;

%% (1) V'*V + L'*L
vars = {'s1'};  dom = [0 1];
ZL = {[0;1]};   ZR = {[0;1]};
params = reshape({sparse(2,2),sparse([2 0;-1 0]),sparse([2 -1;0 0])},[3 1 1]);
Pt = copvar({sopvar(params,struct('out',{vars},'in',{vars}),ZL,ZR,struct('out',dom,'in',dom),[1 1])});
prog = lpiprogram_sop(Pt);
[prog,Pop,inf1] = lpi_ineq_sop(prog,Pt);
prog = lpisolve(prog,sopts);
rel_b = cx_resid(prog);
if rel_b>1e-6
    error('test_lpi_ineq_sop: (1) V''V + L''L is not certified (rel_b %.2e).',rel_b);
end
fprintf('  passed: (1) V''V + L''L certified, rel_b %.1e, lift %s weight %s, terms %s offsets %s, ndec %d\n',...
        rel_b,mat2str(inf1.degrees.D),mat2str(inf1.degrees.w),mat2str(inf1.terms.codes),mat2str(inf1.terms.offsets),numel(Pop.Zd));
npass = npass+1;

%% (2) the 1-D heat LPI through two inequalities
pie = heatNd_pie(1,0);  vars = pie.vars;    dom = pie.dom;    ep = 0.1;
res = zeros(1,2);   kap = [9, 11];
for i = 1:2
    prog = lpiprogram_sop(vars,dom);
    [prog,P] = lpivar_cdopvar(prog,1,{vars},dom,1);
    PT = P'*pie.T;  PA = P'*pie.A0;  TT = pie.T'*pie.T;
    E1 = PT + PT' - 2*ep^2*TT;                  % 2 (P*T - ep^2 T*T), self-adjoint form
    prog = lpi_eq_sop(prog,PT-PT','symmetric');
    [prog,~,iR] = lpi_ineq_sop(prog,E1);
    X = -(PA + PA' + kap(i)*(PT + PT'));        % -(P*A + A*P + 2 kappa P*T) >= 0
    [prog,~,iQ] = lpi_ineq_sop(prog,X);
    prog = lpisolve(prog,sopts);
    res(i) = cx_resid(prog);
    pinf(i) = prog.solinfo.info.pinf; %#ok<AGROW>
    if i==1
        fprintf('  (2) degrees: R lift %s weight %s; Q lift %s weight %s; terms %s offsets %s\n',...
                mat2str(iR.degrees.D),mat2str(iR.degrees.w),mat2str(iQ.degrees.D),mat2str(iQ.degrees.w),...
                mat2str(iQ.terms.codes),mat2str(iQ.terms.offsets));
    end
end
if ~(res(1)<=1e-6 && pinf(1)==0)
    error('test_lpi_ineq_sop: (2) kappa = 9 should be feasible (rel_b %.2e, pinf %d).',res(1),pinf(1));
end
if pinf(2)==0 && res(2)<=1e-6
    error('test_lpi_ineq_sop: (2) kappa = 11 > pi^2 was accepted (rel_b %.2e).',res(2));
end
fprintf('  passed: (2) heat LPI: kappa 9 feasible (rel_b %.1e), kappa 11 rejected (pinf %d, rel_b %.1e)\n',res(1),pinf(2),res(2));
npass = npass+1;

%% (3) 2-D: T'*T >= 0
pie2 = heatNd_pie(2,0);
TT2 = pie2.T'*pie2.T;
prog = lpiprogram_sop(pie2.vars,pie2.dom);
[prog,Pop2,inf3] = lpi_ineq_sop(prog,TT2);
prog = lpisolve(prog,sopts);
rel_b = cx_resid(prog);
if rel_b>1e-6
    error('test_lpi_ineq_sop: (3) 2-D T''T is not certified (rel_b %.2e).',rel_b);
end
fprintf('  passed: (3) 2-D T''T certified, rel_b %.1e, lift %s weight %s, terms %s, ndec %d\n',...
        rel_b,mat2str(inf3.degrees.D),mat2str(inf3.degrees.w),mat2str(inf3.terms.codes),numel(Pop2.Zd));
npass = npass+1;

%% (4) options
prog = lpiprogram_sop(Pt);
[~,Pa] = lpi_ineq_sop(prog,Pt,struct('psatz','none'));
[~,Pb] = lpi_ineq_sop(prog,Pt,struct('psatz','faces'));
[~,Pc,ic] = lpi_ineq_sop(prog,Pt,struct('deg',struct('int',1,'mult',1)));
[~,Pd] = lpi_ineq_sop(prog,Pt,struct('prune',false));
[~,Pe] = lpi_ineq_sop(prog,Pt,struct('dD',0,'dw',0));
if ~(numel(Pa.Zd)<numel(Pop.Zd) && numel(Pb.Zd)>numel(Pop.Zd) && numel(Pe.Zd)<numel(Pop.Zd) && numel(Pd.Zd)>=numel(Pop.Zd))
    error('test_lpi_ineq_sop: (4) the decision counts of the option variants are not ordered as expected.');
end
if ~isfield(ic.degrees,'given')
    error('test_lpi_ineq_sop: (4) a given ''deg'' should bypass the reader.');
end
fprintf('  passed: (4) options: none %d < default %d < faces %d; given deg %d; prune off %d; dD = dw = 0 %d\n',...
        numel(Pa.Zd),numel(Pop.Zd),numel(Pb.Zd),numel(Pc.Zd),numel(Pd.Zd),numel(Pe.Zd));
npass = npass+1;

%% (5) legacy dispatch
pvar s th;
prog1 = lpiprogram(s,th,[0 1]);
[prog1,Pl] = poslpivar(prog1,[1;1],1);
p1 = lpi_ineq(prog1,Pl);
p2 = lpi_ineq_sop(prog1,Pl);
if ~isequal(numel(p1.decvartable),numel(p2.decvartable)) || p1.expr.num~=p2.expr.num
    error('test_lpi_ineq_sop: (5) the legacy dispatch differs from lpi_ineq.');
end
fprintf('  passed: (5) legacy dopvar: lpi_ineq and lpi_ineq_sop agree (%d decision variables, %d expressions)\n',...
        numel(p2.decvartable),p2.expr.num);
npass = npass+1;

%% (6) opts.like: the support of another operator, laid out over P's spaces
Pq = Pt'*Pt;                                % a larger support than Pt's
[~,~,iq] = get_lift_degs(Pq);
[~,Pl,il] = lpi_ineq_sop(prog,Pt,struct('like',Pq));
if ~(isequal(il.degrees.D,iq.D) && isequal(il.degrees.w,iq.w) && numel(Pl.Zd)>numel(Pop.Zd))
    error('test_lpi_ineq_sop: (6) ''like'' should size P from the support of the given operator.');
end
try
    lpi_ineq_sop(prog,Pt,struct('like',1));     e6 = '';
catch ME,   e6 = ME.message;
end
if isempty(e6)
    error('test_lpi_ineq_sop: (6) a non-container ''like'' should be refused.');
end
fprintf('  passed: (6) like: lift %s weight %s from the given operator (own %s %s), ndec %d > %d; non-container refused\n',...
        mat2str(il.degrees.D),mat2str(il.degrees.w),mat2str(inf1.degrees.D),mat2str(inf1.degrees.w),numel(Pl.Zd),numel(Pop.Zd));
npass = npass+1;

fprintf('\ntest_lpi_ineq_sop passed (%d checks).\n',npass);
end
