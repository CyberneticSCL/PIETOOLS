% diag_deg.m -- are redx/redy actually producing a different spec from uniform,
% and does poslpivar_2d even USE dx/dy for this operator?
cuadmm_path;
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+2*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
PIE = initialize(PIE);
Top = PIE.T;
fprintf('DG Top.dim(:,1) = [%s]   <- which components EXIST\n', num2str(Top.dim(:,1)'));

st = lpisettings('heavy');  s2d = st.settings_2d;
dx=s2d.LF_deg.dx; dy=s2d.LF_deg.dy; d2=s2d.LF_deg.d2;
eqd.dx = {1+dx{1}; 1+dx{2}; 1+dx{3}};
eqd.dy = {1+dy{1}, 1+dy{2}, 1+dy{3}};
eqd.d2 = cellfun(@(c) 1+c, d2, 'UniformOutput', false);

fprintf('DG eqd.dx{1}=%s  eqd.dy{1}=%s\n', mat2str(eqd.dx{1}), mat2str(eqd.dy{1}));
fprintf('DG eqd.d2{1,1} =\n'); disp(eqd.d2{1,1});
A = redx(eqd,1); B = redy(eqd,1); C = bump(eqd,-1);
fprintf('DG redx.d2{1,1} =\n'); disp(A.d2{1,1});
fprintf('DG redy.d2{1,1} =\n'); disp(B.d2{1,1});
fprintf('DG unif.d2{1,1} =\n'); disp(C.d2{1,1});
fprintf('DG dx after redx: %s   dy after redx: %s\n', mat2str(A.dx{1}), mat2str(A.dy{1}));
fprintf('DG dx after unif: %s   dy after unif: %s\n', mat2str(C.dx{1}), mat2str(C.dy{1}));
fprintf('DG d2{2,2} eqd:\n'); disp(eqd.d2{2,2});
fprintf('DG d2{2,2} redx:\n'); disp(A.d2{2,2});
fprintf('DG d2{2,2} unif:\n'); disp(C.d2{2,2});

% now: do the resulting Gram blocks actually differ?
Aop = PIE.A;
np = Top.dim(:,1);
Iop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop = clean_opvar((Aop'*(Iop*Top))' + Aop'*(Iop*Top),1e-12);
eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
fprintf('DG eq_opts.exclude = %s\n', mat2str(double(eq_opts.exclude)));
for nm = {'eqd','redx','redy','unif'}
    switch nm{1}
        case 'eqd',  dd = eqd; case 'redx', dd = A;
        case 'redy', dd = B;   case 'unif', dd = C;
    end
    p2 = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
    o = eq_opts; o.psatz = 3;
    [p2,Qe] = poslpivar_2d(p2,Qop.dim,dd,o);
    S = cuadmm_private('sdpshape',p2);
    fprintf('DG psatz3 with %-5s -> block %s\n', nm{1}, mat2str(S.Ks));
end
fprintf('DGDONE\n');

function d = redx(d,k)
d.dx = cellfun(@(c) max(c-k,0), d.dx, 'UniformOutput', false);
d.d2 = cellfun(@(M) rowcut(M,k), d.d2, 'UniformOutput', false);
end
function d = redy(d,k)
d.dy = cellfun(@(c) max(c-k,0), d.dy, 'UniformOutput', false);
d.d2 = cellfun(@(M) colcut(M,k), d.d2, 'UniformOutput', false);
end
function M = rowcut(M,k), if size(M,1)>1, M(2:end,:) = max(M(2:end,:)-k,0); end, end
function M = colcut(M,k), if size(M,2)>1, M(:,2:end) = max(M(:,2:end)-k,0); end, end
function d = bump(d,k)
d.dx = cellfun(@(c) max(c+k,0), d.dx, 'UniformOutput', false);
d.dy = cellfun(@(c) max(c+k,0), d.dy, 'UniformOutput', false);
d.d2 = cellfun(@(c) max(c+k,0), d.d2, 'UniformOutput', false);
end
