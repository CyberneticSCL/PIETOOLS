% psatz1d.m -- the same generator question in 1-D, with P pinned to I.
%
% On [0,1] the choices are {s(1-s)} (what PIETOOLS ships), {s, 1-s} as two
% SEPARATE multipliers, or all three.  2-D showed that separate LINEAR
% generators beat the product on BOTH reach and cost; the question here is
% whether the same holds in 1-D, and specifically whether the right generators
% buy the reach that currently costs a degree bump (Dup 1 -> 2, m 32 -> 53).
% lam* = pi^2 = 9.8696.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('poslpivar');
fprintf('Q1 poslpivar resolves to %s\n', which('poslpivar'));

LAMSTAR = pi^2;
FR = [0.10 0.50 0.90 0.99 0.999 0.9999];

% codes: 1 = s(1-s) product (stock), 2 = (s-a), 3 = (b-s)
V = { 'Dup1 none'              , 1, []
      'Dup1 product {s(1-s)}'  , 1, 1
      'Dup1 linear {s,1-s}'    , 1, [2 3]
      'Dup1 all three'         , 1, [2 3 1]
      'Dup2 product (current)' , 2, 1
      'Dup2 linear {s,1-s}'    , 2, [2 3] };

fprintf('Q1 variant|m|nblk|Ns|frac|lam|cert|rel_b|feasratio|t\n');
for v = 1:size(V,1)
    lab = V{v,1}; Dup = V{v,2}; codes = V{v,3}; best = 0; Sm = [];
    for f = FR
        lam = f*LAMSTAR;
        try
            t0=tic; prog = build1(lam,Dup,codes);
            so.solver='sedumi'; so.params.fid=0;
            evalc('sol = lpisolve(prog,so);'); tw=toc(t0);
            S = cuadmm_private('sdpshape',sol); if isempty(Sm), Sm=S; end
            I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6) || norm(xv)<=1e-12;
            cert = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
            fprintf('Q1 %s|%d|%d|[%s]|%.4f|%.4f|%d|%.3e|%+.4f|%.2f\n', ...
                lab,S.m,S.nblk,strtrim(num2str(S.Ks)),f,lam,cert,rb,I.feasratio,tw);
            if cert, best=f; else, break; end
        catch ME
            fprintf('Q1 %s|ERR|%.4f|%s\n',lab,f,strrep(ME.message,newline,' ')); break
        end
    end
    if isempty(Sm), Sm.m=-1; Sm.nblk=-1; end
    fprintf('Q1REACH %s|reach=%.4f|m=%d|nblk=%d\n',lab,best,Sm.m,Sm.nblk);
end
fprintf('Q1DONE\n');

cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function prog = build1(lam,Dup,codes)
pvar s t
x = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+lam*x; subs(x,s,0)==0; subs(x,s,1)==0]);
PIE = initialize(PIE);
Top = PIE.T; Aop = PIE.A;
st = lpisettings('light'); st.eppos2 = 1e-2; st.eppos = 1e-2;
n1=1;n2=1;n3=1;n4=n2+n3;
dd2 = {n1+Dup,   [n2+Dup-1, n3+Dup,   n4+Dup  ], [n2+Dup-1, n3+Dup,   n4+Dup  ]};
dd3 = {n1+Dup-1, [n2+Dup-2, n3+Dup-1, n4+Dup-1], [n2+Dup-2, n3+Dup-1, n4+Dup-1]};

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
opvar Iop; Iop.I = PIE.dom; Iop.var1 = PIE.vars(1,1); Iop.var2 = PIE.vars(1,2);
Iop.dim = Top.dim;
if Top.dim(1,1)>0, Iop.P = eye(Top.dim(1,1)); end
if Top.dim(2,1)>0, Iop.R.R0 = eye(Top.dim(2,1)); end
Dop = clean_opvar(Aop'*(Iop*Top) + (Iop*Top)'*Aop, 1e-12);

o0 = st.options2; o0.psatz = 0;
[prog,Deop] = poslpivar(prog,Dop.dim,dd2,o0);
for c = codes
    o = st.options3; o.psatz = c;
    [prog,De2] = poslpivar(prog,Dop.dim,dd3,o);
    Deop = Deop + De2;
end
prog = lpi_eq(prog,Dop+Deop,'symmetric');
end
