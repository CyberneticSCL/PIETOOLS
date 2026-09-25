% slot101.m -- verify the corrected 1-D psatz reduction MYSELF, on both axes.
%
% SOURCE FACT (verified by reading poslpivar.m:342-343): the substitutions
%   Z2etath = subs(Z2sth,var1,sss);  Z2etas = subs(Z2etath,var2,var1);
% send the spec's FIRST slot to the INTEGRATION variable and the SECOND to the
% SURVIVING kernel variable -- the opposite of the field names.  The psatz
% multiplier is evaluated at the quadrature variable, so the reduction belongs
% in slot 1 and the joint slot, NOT slot 2.  Stock uses a uniform [1 1 1].
%
% TEST 1 (support): does dd2-[1 0 1] reproduce the base block's assembled
%   monomial support where stock dd2-[1 1 1] loses some?
% TEST 2 (reach): nested specs give nested cones, so the larger cone must have
%   reach >= stock.  Does it actually BUY anything?  The workflow only checked
%   support; reach is the thing that matters and was never measured.
cuadmm_path;
pvar s t
LAMSTAR = pi^2;

% ---------- TEST 1: assembled monomial support ----------
fprintf('SL --- support at Dup=1 and Dup=2\n');
for Dup = [1 2]
    [dd2,dd3u,dd3r] = specs(Dup);
    base = supp(dd2 ,0);
    stok = supp(dd3u,1);
    rule = supp(dd3r,1);
    fprintf('SL Dup=%d base=%d | stock(-[1 1 1])=%d missing=%d extra=%d | rule(-[1 0 1])=%d missing=%d extra=%d\n', ...
        Dup, size(base,1), ...
        size(stok,1), nmiss(base,stok), nmiss(stok,base), ...
        size(rule,1), nmiss(base,rule), nmiss(rule,base));
end

% ---------- TEST 2: reach ----------
fprintf('SL --- reach (P=I, SeDuMi)\n');
FR = [0.50 0.75 0.90 0.99 0.999 0.9999];
for Dup = [1 2]
  for mode = {'stock','rule'}
    best = 0; mm = -1; nb = -1;
    for f = FR
        try
            [prog,S] = build(Dup,mode{1},f*LAMSTAR);
            so.solver='sedumi'; so.params.fid=0;
            evalc('sol = lpisolve(prog,so);');
            if mm<0, mm=S.m; nb=S.nblk; ks=S.Ks; end
            I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
            cert = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
            if cert, best=f; else, break; end
        catch ME
            fprintf('SL Dup=%d %s ERR %s\n',Dup,mode{1},strrep(ME.message,newline,' ')); break
        end
    end
    fprintf('SL Dup=%d %-5s reach=%.4f m=%d nblk=%d Ns=[%s]\n',Dup,mode{1},best,mm,nb,strtrim(num2str(ks)));
  end
end
fprintf('SLDONE\n');

function [dd2,dd3u,dd3r] = specs(Dup)
n1=1;n2=1;n3=1;n4=2;
dd2  = {n1+Dup, [n2+Dup-1, n3+Dup, n4+Dup], [n2+Dup-1, n3+Dup, n4+Dup]};
% stock: uniform -1 on every entry
dd3u = {n1+Dup-1, [n2+Dup-2, n3+Dup-1, n4+Dup-1], [n2+Dup-2, n3+Dup-1, n4+Dup-1]};
% rule: -floor(e/2)*[1 0 1] with e=2 for s(1-s), i.e. -1 on the QUADRATURE slot
% and the JOINT slot, and 0 on the SURVIVING slot.
dd3r = {n1+Dup-1, [n2+Dup-2, n3+Dup,   n4+Dup-1], [n2+Dup-2, n3+Dup,   n4+Dup-1]};
end

function M = supp(dd,ps)
% assembled monomial support of the R1 kernel of a poslpivar block
pvar s t
x = pde_var('state',1,s,[0,1]);
PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+2*x; subs(x,s,0)==0; subs(x,s,1)==0]));
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
o.psatz = ps; o.exclude = [0 0 0 0]; o.sep = 0;
[~,P] = poslpivar(prog,PIE.T.dim,dd,o);
R1 = P.R.R1;
if isempty(R1) || isempty(R1.degmat), M = zeros(0,2); return; end
M = unique(full(R1.degmat),'rows');
end

function n = nmiss(A,B)
% rows of A not present in B
n = 0;
for i=1:size(A,1)
    if ~any(all(B==A(i,:),2)), n = n+1; end
end
end

function [prog,S] = build(Dup,mode,lam)
pvar s t
x = pde_var('state',1,s,[0,1]);
PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+lam*x; subs(x,s,0)==0; subs(x,s,1)==0]));
Top = PIE.T; Aop = PIE.A;
[dd2,dd3u,dd3r] = specs(Dup);
if strcmp(mode,'stock'), dd3 = dd3u; else, dd3 = dd3r; end
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
opvar Iop; Iop.I = PIE.dom; Iop.var1 = PIE.vars(1,1); Iop.var2 = PIE.vars(1,2);
Iop.dim = Top.dim;
if Top.dim(1,1)>0, Iop.P = eye(Top.dim(1,1)); end
if Top.dim(2,1)>0, Iop.R.R0 = eye(Top.dim(2,1)); end
Dop = clean_opvar(Aop'*(Iop*Top) + (Iop*Top)'*Aop, 1e-12);
o0.psatz=0; o0.exclude=[0 0 0 0]; o0.sep=0;
o1.psatz=1; o1.exclude=[0 0 0 0]; o1.sep=0;
[prog,D1] = poslpivar(prog,Dop.dim,dd2,o0);
[prog,D2] = poslpivar(prog,Dop.dim,dd3,o1);
prog = lpi_eq(prog,Dop+D1+D2,'symmetric');
S = cuadmm_private('sdpshape',prog);
end
