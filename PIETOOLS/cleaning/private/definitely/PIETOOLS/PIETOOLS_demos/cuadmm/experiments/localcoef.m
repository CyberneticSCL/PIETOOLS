% localcoef.m -- LOCALLY positive coefficients, which is where a
% Positivstellensatz certificate should actually earn its keep.
%
% Maintainer's point: variable coefficients need psatz MORE than homogeneous
% PDEs, because the coefficients may only be locally positive.  My earlier
% variable-coefficient plants used a(s) = 1+s, positive everywhere on [0,1]
% (range 1..2) -- globally positive, so they never stressed locality at all.
%
% Plants here make the coefficient's positivity GENUINELY LOCAL:
%   L1  destabilising reaction localised near s1=0:  lam*(1-2*s1)*x
%       (the reaction is +lam at s1=0 and -lam at s1=1, so the instability
%        lives on half the domain and the certificate must be local)
%   L2  near-degenerate diffusion: a(s1) = eps + s1(1-s1), which VANISHES at
%       both s1 faces as eps->0.  Written in divergence form
%       d/ds1(a dx/ds1) = a*x_s1s1 + (1-2*s1)*x_s1.
%   L3  reaction positive only on a central patch: lam*(1-4*(s1-0.5)^2)*x
%   L0  constant control at the same lambda ladder, for reference.
%
% PREDICTION UNDER TEST (the maintainer's): the psatz-vs-no-psatz GAP should
% WIDEN relative to the homogeneous case.  Within-plant comparison; lam* is
% not analytic for any of these.
cuadmm_path;
LAM  = [1 2 5 10 15 19];
SETS = {'none',[] ; 'product',1 ; 'linear4',[3 4 5 6]};
P = {'L0 constant',0; 'L1 local destab reaction',1; 'L2 degenerate diffusion',2; 'L3 central patch',3};

fprintf('LC plant|set|m|nblk|lam|cert|rel_b|feasratio|numerr\n');
for ip = 1:size(P,1)
  for is = 1:size(SETS,1)
    nm=SETS{is,1}; codes=SETS{is,2}; best=-inf; mm=-1; nb=-1;
    for lam = LAM
      try
        prog = b10(P{ip,2},lam,codes);
        so.solver='mosek'; evalc('sol = lpisolve(prog,so);');
        S=cuadmm_private('sdpshape',sol); if mm<0, mm=S.m; nb=S.nblk; end
        I=sol.solinfo.info; xv=sol.solinfo.RRx(:);
        Atf=[];bf=[];
        for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
        rb=norm(full(Atf'*xv-bf))/norm(full(bf));
        triv=(abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
        cert=(I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
        fprintf('LC %s|%s|%d|%d|%g|%d|%.3e|%+.4f|%d\n', ...
                P{ip,1},nm,S.m,S.nblk,lam,cert,rb,I.feasratio,I.numerr);
        if cert, best=lam; else, break; end
      catch ME
        fprintf('LC %s|%s|ERR|%g|%s\n',P{ip,1},nm,lam,strrep(ME.message,newline,' ')); break
      end
    end
    if isinf(best), bs='none'; else, bs=sprintf('%g',best); end
    fprintf('LCREACH %s|%s|lam_max=%s|m=%d|nblk=%d\n',P{ip,1},nm,bs,mm,nb);
  end
end
fprintf('LCDONE\n');

function prog = b10(kind,lam,codes)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0 1;0 1]);
switch kind
  case 0, rhs = diff(x,s1,2)+diff(x,s2,2)+lam*x;
  case 1, rhs = diff(x,s1,2)+diff(x,s2,2)+lam*(1-2*s1)*x;
  case 2, a = 0.05+s1*(1-s1);
          rhs = a*diff(x,s1,2)+(1-2*s1)*diff(x,s1,1)+diff(x,s2,2)+lam*x;
  case 3, rhs = diff(x,s1,2)+diff(x,s2,2)+lam*(1-4*(s1-0.5)^2)*x;
end
PIE = initialize(convert([diff(x,t,1)==rhs;
      subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0]));
Top=PIE.T; Aop=PIE.A;
st=lpisettings('heavy'); st.settings_2d.eppos=1e-2*[1;1;1;1];
s2d=st.settings_2d;
dx=s2d.LF_deg.dx; dy=s2d.LF_deg.dy; d2=s2d.LF_deg.d2;
eqd.dx={1+dx{1};1+dx{2};1+dx{3}};
eqd.dy={1+dy{1},1+dy{2},1+dy{3}};
eqd.d2=cellfun(@(c) 1+c, d2,'UniformOutput',false);
prog=lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
np=Top.dim(:,1);
Iop=opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop=clean_opvar((Aop'*(Iop*Top))'+Aop'*(Iop*Top),1e-12);
eq_opts=get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop]=poslpivar_2d(prog,Qop.dim,eqd,eq_opts);
for p=codes
    o=eq_opts; o.psatz=p;
    [prog,Qe]=poslpivar_2d(prog,Qop.dim,eqd,o);
    Qeop=Qeop+Qe;
end
prog=lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end
