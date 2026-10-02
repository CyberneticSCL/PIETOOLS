function R = test_copquadvar_faces()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_COPQUADVAR_FACES() asserting tests of the face codes of
% copquadvar / poscopvar (options.psatz = 2d+1 / 2d+2, MMP 09/27/2026), the
% Psatz terms of HEATND_LPI 'linear' and HEATND_POINCARE 'faces':
%
%  (i)   value checks, N = 1, 2, 3: psatz 2 (poslpivar_2d's ball), 2N+3,
%        1.5, [3 4], and type 'sym' with a face code error; 0, 1 and every
%        code 3..2N+2 are accepted; int32/uint8 codes build exactly the
%        operator of the double code;
%  (ii)  semantics, from the definition: at degree 0 each basis operator is
%        (Z_alpha x)(th) = int I_alpha(th - s) x(s) ds (alpha_k = 1 delta,
%        2 s_k <= th_k, 3 s_k >= th_k), so for a random Q >= 0 put into the
%        declared variables (by Qcell),
%            <x, Pop x>  =  int g(th) (Z x)(th)' Q (Z x)(th) dth,
%        the left side by HEATND_APPLY of the substituted operator, the
%        right side by independent quadrature, with g built HERE from the
%        header's definition: code 2d+1 (th_d-a_d)/L_d, 2d+2 (b_d-th_d)/L_d
%        (d over the sorted registry), 1 the product, 0 none. N = 1, 2, 3 on
%        a non-unit box, every code. A weight at s instead of th, the wrong
%        direction, or lower/upper swapped fails (ii);
%  (iii) control: the right side with the OTHER face of the same direction
%        must NOT match (so (ii) can see a lower/upper swap);
%  (iv)  the Galerkin matrix of Pop on the tensor monomials of degree <= 2
%        per variable is symmetric PSD (N = 1, 2).
%
% Initial coding MMP, 09/27/2026 (from heatNd_test_posw, retired with the
%   copy heatNd_posw it tested; (ii) is its check (ii), extended to N = 3
%   and to every code).
% MMP, 10/01/2026: The program is lpiprogram for every N; it no longer refuses
%   N > 2, so the hand-built copy for N > 2 is commented out. Same program
%   (lpiprogram builds exactly what the copy built).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

heatNd_path();
rng(4242);
w0 = warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
R = struct('check',{},'err',{},'tol',{});
Dnu = [0 1; -1 2; 0.5 3];
% % % (i) value checks.
for N = 1:3
    vars = arrayfun(@(i) sprintf('s%d',i),1:N,'UniformOutput',false);
    dom = Dnu(1:N,:);
    for bad = {2, 2*N+3, 1.5, [3 4]}
        ok = errs(@() poscopvar(mkprog(vars,dom),1,vars,dom,0,struct('psatz',bad{1})));
        R(end+1) = struct('check',sprintf('(i) N=%d psatz=%s rejected',N,mat2str(bad{1})),'err',double(~ok),'tol',0); %#ok<AGROW>
    end
    ok = errs(@() copquadvar(mkprog(vars,dom),1,vars,dom,0,struct('psatz',2*N+2,'type','sym')));
    R(end+1) = struct('check',sprintf('(i) N=%d psatz=%d, type sym rejected',N,2*N+2),'err',double(~ok),'tol',0); %#ok<AGROW>
    for c = [0 1 3:2*N+2]
        ok = ~errs(@() poscopvar(mkprog(vars,dom),1,vars,dom,0,struct('psatz',c)));
        R(end+1) = struct('check',sprintf('(i) N=%d psatz=%d accepted',N,c),'err',double(~ok),'tol',0); %#ok<AGROW>
    end
    % Integer classes pass the value check, so they must build the operator
    % of the double code (int32(4) once built the upper face of s2 at N=2).
    for c = 3:2*N+2
        [~,P0] = poscopvar(mkprog(vars,dom),1,vars,dom,0,struct('psatz',c));
        for cls = {@int32,@uint8}
            [~,P1] = poscopvar(mkprog(vars,dom),1,vars,dom,0,struct('psatz',cls{1}(c)));
            R(end+1) = struct('check',sprintf('(i) N=%d psatz=%s(%d) equals double',N,func2str(cls{1}),c), ...
                              'err',double(~isequal(P0,P1)),'tol',0); %#ok<AGROW>
        end
    end
end
% % % (ii)-(iv) semantics at degree 0.
for N = 1:3
    vars = arrayfun(@(i) sprintf('s%d',i),1:N,'UniformOutput',false);
    dom = Dnu(1:N,:);
    nq = 8;     if N==3,    nq = 6;     end     % exact: integrands of degree <= 2nq-1 per direction
    [Xq,wq] = grid(dom,nq);
    for code = [0 1 3:2*N+2]
        [g,gw] = weight(code,dom);
        [~,Pop,Qcell,bl] = poscopvar(mkprog(vars,dom),1,vars,dom,0,struct('psatz',code));
        nb = size(bl,1);
        G = randn(nb);  Q = G*G';                   % random PSD
        W = subsd(Pop,Qcell,Q);
        Pw = randi([0 2],4,N);  cf = randn(4,1);    % test function, degree <= 2 per variable
        x = @(X) sum(cf'.*prod(reshape(X,[],1,N).^reshape(Pw,1,[],N),3),2);
        lhs = sum(wq.*x(Xq).*heatNd_apply(W,x,Xq,nq));
        Zx = zeros(size(Xq,1),nb);                  % (Z_alpha x)(th) at th = Xq
        for c = 1:nb,   Zx(:,c) = basis_apply(bl(c,2:end),x,Xq,dom,nq);  end
        q = sum((Zx*Q).*Zx,2);
        rhs = sum(wq.*g(Xq).*q);
        R(end+1) = struct('check',sprintf('(ii) N=%d psatz=%d: <x,Wx> = int g (Zx)''Q(Zx)',N,code), ...
                          'err',abs(lhs-rhs)/max(abs(rhs),eps),'tol',1e-10); %#ok<AGROW>
        if ~isempty(gw)
            rw = sum(wq.*gw(Xq).*q);
            R(end+1) = struct('check',sprintf('(iii) N=%d psatz=%d: other face does not match',N,code), ...
                              'err',double(abs(lhs-rw)/max(abs(rw),eps)<1e-3),'tol',0); %#ok<AGROW>
        end
        if N<=2
            E = dec2base(0:3^N-1,3,N)-'0';      nphi = size(E,1);
            F = zeros(size(Xq,1),nphi);     WF = F;
            for j = 1:nphi
                phi = @(X) prod(X.^E(j,:),2);
                F(:,j) = phi(Xq);   WF(:,j) = heatNd_apply(W,phi,Xq,nq);
            end
            Gm = F'*(wq.*WF);   ev = eig((Gm+Gm')/2);
            e4 = max(norm(Gm-Gm','fro')/norm(Gm,'fro'), max(0,-min(ev)/max(abs(ev))));
            R(end+1) = struct('check',sprintf('(iv) N=%d psatz=%d: Galerkin %dx%d symmetric PSD',N,code,nphi,nphi), ...
                              'err',e4,'tol',1e-10);    %#ok<AGROW>
        end
    end
end
warning(w0);
ok = true;
for k = 1:numel(R)
    pass = R(k).err<=R(k).tol;     ok = ok && pass;
    fprintf('%-62s %10.2e  %s\n',R(k).check,R(k).err,char(string(pass)));
end
assert(ok,'test_copquadvar_faces: a check failed.');
fprintf('test_copquadvar_faces: all %d checks pass.\n',numel(R));
end


function tf = errs(f)
% True if F() errors.
tf = false;
try,    f();    catch,  tf = true;  end
end

function [g,gw] = weight(code,dom)
% The weight of CODE from copquadvar's header: 0 none, 1 the product,
% 2d+1 (th_d-a_d)/L_d, 2d+2 (b_d-th_d)/L_d. GW: the other face of d.
gw = [];
switch code
    case 0,     g = @(T) ones(size(T,1),1);                         return
    case 1,     g = @(T) prod((T-dom(:,1)').*(dom(:,2)'-T),2);      return
end
d = floor((code-1)/2);  L = dom(d,2)-dom(d,1);
lo = @(T) (T(:,d)-dom(d,1))/L;      hi = @(T) (dom(d,2)-T(:,d))/L;
if mod(code,2),     g = lo;     gw = hi;
else,               g = hi;     gw = lo;
end
end

function prog = mkprog(vars,dom)
% lpiprogram refuses N > 2; for N = 3 what it would return (as HEATND_LPI).
% (10/01/2026: it no longer refuses; one call for every N.)                 % MMP, 10/01/2026
% if numel(vars)<=2                                                         % MMP, 10/01/2026 (was)
%     prog = lpiprogram(polynomial(vars(:)),[],dom);  return                % MMP, 10/01/2026 (was)
% end                                                                       % MMP, 10/01/2026 (was)
% prog = sosprogram(polynomial([]),dpvar(zeros(0,1)));                      % MMP, 10/01/2026 (was)
% prog.vartable = [prog.vartable; polynomial(vars(:)); polynomial(strcat(vars(:),'_dum'))]; % MMP, 10/01/2026 (was)
% prog.dom = dom;                                                           % MMP, 10/01/2026 (was)
prog = lpiprogram(polynomial(vars(:)),[],dom);                              % MMP, 10/01/2026
end

function y = basis_apply(alpha,x,Th,dom,nq)
% (Z_alpha x)(th) = int I_alpha(th - s) x(s) ds, degree-0 basis (kernel 1).
[xg,wg] = gl01(nq);     N = numel(alpha);   y = zeros(size(Th,1),1);
for k = 1:size(Th,1)
    th = Th(k,:);   pts = cell(1,N);    wts = cell(1,N);
    for d = 1:N
        switch alpha(d)
            case 1,     pts{d} = th(d);     wts{d} = 1;
            case 2,     pts{d} = dom(d,1)+(th(d)-dom(d,1))*xg;  wts{d} = (th(d)-dom(d,1))*wg;
            case 3,     pts{d} = th(d)+(dom(d,2)-th(d))*xg;     wts{d} = (dom(d,2)-th(d))*wg;
            otherwise,  error('basis_apply:alpha','Unexpected multi-index %d.',alpha(d));
        end
    end
    [P,Wt] = tens(pts,wts);
    y(k) = sum(Wt.*x(P));
end
end

function [X,w] = grid(dom,n)
[xg,wg] = gl01(n);  N = size(dom,1);    pts = cell(1,N);    wts = cell(1,N);
for d = 1:N,    L = dom(d,2)-dom(d,1);  pts{d} = dom(d,1)+L*xg;  wts{d} = L*wg;  end
[X,w] = tens(pts,wts);
end

function [X,w] = tens(pts,wts)
% Tensor grid of per-direction nodes/weights; columns in direction order.
N = numel(pts);     X = zeros(1,0);     w = 1;
for d = 1:N
    n = numel(pts{d});  m = size(X,1);
    X = [repmat(X,n,1), repelem(pts{d}(:),m,1)];
    w = repmat(w,n,1).*repelem(wts{d}(:),m,1);
end
end

function Y = subsd(X,Qcell,Q)
% Fixed operator X(d) of a 1x1 'cdopvar' with Qcell{c,e} := Q(c,e), other
% decision variables 0: params{g} = unvec(A_g + B_g'*d) (sdopvar.m).
names = {};     vals = [];
for c = 1:size(Qcell,1)
    for e = 1:size(Qcell,2)
        nm = cellstr(string(Qcell{c,e}));  names = [names; nm(:)];          %#ok<AGROW>
        vals = [vals; Q(c,e)*ones(numel(nm),1)];                            %#ok<AGROW>
    end
end
B = X.C{1};     zd = B.Zd;  if ~iscellstr(zd),  zd = cellstr(string(zd));  end
[tf,loc] = ismember(zd(:),names(:));
d = zeros(numel(zd),1);     d(tf) = vals(loc(tf));
nL = prod(cellfun(@numel,B.ZL));    nR = prod(cellfun(@numel,B.ZR));
m = B.dims(1)*nL;   n = B.dims(2)*nR;
prm = cell(size(B.params.A));
for g = 1:numel(prm)
    a = B.params.A{g};  Bg = B.params.B{g};
    if isempty(a),  a = sparse(m*n,1);  end
    if ~isempty(Bg),    a = a + Bg'*d;  end
    prm{g} = reshape(full(a),m,n);
end
Y = copvar({sopvar(prm,B.vars,B.ZL,B.ZR,B.dom,B.dims)});
end

function [x,w] = gl01(n)
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     w = 2*V(1,ix)'.^2;
x = (x+1)/2;    w = w/2;
end
