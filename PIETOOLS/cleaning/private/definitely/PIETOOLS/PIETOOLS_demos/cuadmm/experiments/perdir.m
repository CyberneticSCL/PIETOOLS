% perdir.m -- THE decisive test: PER-DIRECTION psatz degree reduction.
%
% Theory (derived, not swept).  poslpivar_2d builds Pop = Z'((Q*g)Z), so the
% x-degree of the result is 2*dx + e_x where e_x = deg_x(g).  Matching a target
% of x-degree Dx gives dx = (Dx - e_x)/2 = dx0 - e_x/2.  For the LINEAR
% generators e = 1, so no integer choice is exactly balanced: dx0 is over by
% one, dx0 - 1 is under by one.  Over-degree forces the top-degree part of
% g*Sigma to vanish; since Sigma is PSD that pins its top diagonal block, and
% hence its whole row/column, to zero -- truncating the effective basis to
% dx0 - 1.  PREDICTION: dx0 and dx0-1 (applied PER DIRECTION) give the SAME
% cone, so the same reach, with dx0-1 cheaper in m and better conditioned.
%
% From poslpivar_2d.m:66-78 the d2{i,j} arrays are laid out
%        0 |      y |     nu |    y*nu
%        x |    x*y |   x*nu |  x*y*nu
%       tt |   tt*y |  tt*nu | tt*y*nu
%     x*tt | x*tt*y |x*tt*nu |x*tt*y*nu
% so ROWS index the x-direction and COLUMNS the y-direction.  A multiplier in
% s1 reduces dx and rows 2:end; a multiplier in s2 reduces dy and cols 2:end.
% Row 1 / column 1 carry no variable of that direction and must NOT move.
% My earlier uniform subtraction hit dy and row 1 as well, which is why it
% failed -- that tested "uniform reduction", not "reduction".
cuadmm_path;
LAMSTAR = 2*pi^2;
FR = [0.90 0.99 0.999 0.9999];

V = { 'A full eq_deg (baseline)'   , 'full'
      'B per-direction -1 (THEORY)', 'perdir'
      'C uniform -1 (confounded)'  , 'uniform' };

fprintf('PDIR variant|m|nblk|Ns|frac|cert|rel_b|feasratio|numerr|t\n');
for v = 1:size(V,1)
    lab = V{v,1}; mode = V{v,2}; best = 0; Sm = [];
    for f = FR
        lam = f*LAMSTAR;
        try
            t0=tic; prog = build5(lam,mode);
            so.solver='mosek'; evalc('sol = lpisolve(prog,so);'); tw=toc(t0);
            S = cuadmm_private('sdpshape',sol); if isempty(Sm), Sm=S; end
            I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
            cert = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
            fprintf('PDIR %s|%d|%d|[%s]|%.4f|%d|%.3e|%+.4f|%d|%.1f\n', ...
                lab,S.m,S.nblk,strtrim(num2str(S.Ks)),f,cert,rb,I.feasratio,I.numerr,tw);
            if cert, best=f; else, break; end
        catch ME
            fprintf('PDIR %s|ERR|%.4f|%s\n',lab,f,strrep(ME.message,newline,' ')); break
        end
    end
    if isempty(Sm), Sm.m=-1; Sm.nblk=-1; end
    fprintf('PDIRREACH %s|reach=%.4f|m=%d|nblk=%d\n',lab,best,Sm.m,Sm.nblk);
end
fprintf('PDIRDONE\n');

function prog = build5(lam,mode)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
PIE = initialize(PIE);
Top = PIE.T; Aop = PIE.A;
st = lpisettings('heavy'); st.settings_2d.eppos = 1e-2*[1;1;1;1];
s2d = st.settings_2d;
dx=s2d.LF_deg.dx; dy=s2d.LF_deg.dy; d2=s2d.LF_deg.d2;
eqd.dx = {1+dx{1}; 1+dx{2}; 1+dx{3}};
eqd.dy = {1+dy{1}, 1+dy{2}, 1+dy{3}};
eqd.d2 = cellfun(@(c) 1+c, d2, 'UniformOutput', false);

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
np = Top.dim(:,1);
Iop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop = clean_opvar((Aop'*(Iop*Top))' + Aop'*(Iop*Top),1e-12);
eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,eqd,eq_opts);

% psatz codes 3,4 act in s1; codes 5,6 act in s2.
for p = [3 4 5 6]
    o = eq_opts; o.psatz = p;
    switch mode
        case 'full',    pd = eqd;
        case 'uniform', pd = bump(eqd,-1);
        case 'perdir'
            if any(p==[3 4]), pd = redx(eqd,1); else, pd = redy(eqd,1); end
    end
    [prog,Qe] = poslpivar_2d(prog,Qop.dim,pd,o);
    Qeop = Qeop + Qe;
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end

function d = redx(d,k)
% x-direction only: dx and ROWS 2:end of every d2 block.  dy untouched.
d.dx = cellfun(@(c) max(c-k,0), d.dx, 'UniformOutput', false);
d.d2 = cellfun(@(M) rowcut(M,k), d.d2, 'UniformOutput', false);
end
function d = redy(d,k)
% y-direction only: dy and COLUMNS 2:end of every d2 block.  dx untouched.
d.dy = cellfun(@(c) max(c-k,0), d.dy, 'UniformOutput', false);
d.d2 = cellfun(@(M) colcut(M,k), d.d2, 'UniformOutput', false);
end
function M = rowcut(M,k)
if size(M,1)>1, M(2:end,:) = max(M(2:end,:)-k,0); end
end
function M = colcut(M,k)
if size(M,2)>1, M(:,2:end) = max(M(:,2:end)-k,0); end
end
function d = bump(d,k)
d.dx = cellfun(@(c) max(c+k,0), d.dx, 'UniformOutput', false);
d.dy = cellfun(@(c) max(c+k,0), d.dy, 'UniformOutput', false);
d.d2 = cellfun(@(c) max(c+k,0), d.d2, 'UniformOutput', false);
end
