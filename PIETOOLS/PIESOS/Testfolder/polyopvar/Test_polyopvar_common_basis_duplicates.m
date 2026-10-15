% Test_polyopvar_common_basis_duplicates.m -- common_basis (hence
% vertcat, plus) for 'polyopvar' objects that are not canonical: several
% terms with the same monomial (equal rows of Z). The operators acting on
% equal monomials must be added, not overwritten.
%
% A non-canonical object F = f + g (f, g with the same Z, built as one
% object with Z=[Z;Z], C=[Cf,Cg]) has the value eval(f)+eval(g). Oracle:
% the value of each operand of common_basis, and of the rows of vertcat,
% evaluated at random polynomial states by eval_polyopvar, equals that
% before the operation.
% Cases: duplicates in the first operand, in the second, in both;
% monomials shared by the operands; functions (pointwise and Kronecker
% output, 1-2 spatial variables) and functionals (intop, degree one).
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;
ncases = 0;
allnames = {'a','b','c'};

%% Functions
for N = 1:2
    for kind = {'pw','kron'}
        rng(100*N+numel(kind{1}));
        [svars,dom] = spatial_space(N);
        svs = pv_cell(N);
        x = {rand_pol(svs),rand_pol(svs),rand_pol(svs)};
        % Zs{k}: monomials of the operands, and number of terms of each
        specs = {   {[1,2],[1,2]},                  {[1,2]}                 % duplicates in A, overlap with B
                    {[1,2]},                        {[1,2],[1,2]}           % duplicates in B
                    {[1,2],[1,2],[2,3]},            {[2,3],[2,3],[1,2]} };  % both, several monomials
        for k=1:size(specs,1)
            A = raw_poly(N,specs{k,1},kind{1});             B = raw_poly(N,specs{k,2},kind{1});
            lbl = sprintf('function N=%d,%s,case %d',N,kind{1},k);
            VA = ev(A,x,allnames);                          VB = ev(B,x,allnames);
            % the oracle: the sum of the values of the terms
            assert_close(VA,sum_terms(A,x,allnames),tol,[lbl,' (oracle A)']);
            [A2,B2] = common_basis(A,B);
            assert_close(ev(A2,x,allnames),VA,tol,[lbl,' (common_basis A)']);
            assert_close(ev(B2,x,allnames),VB,tol,[lbl,' (common_basis B)']);
            ncases = ncases+3;
            % vertcat, plus
            V = [A;B];
            assert_close(ev(V,x,allnames),[VA;VB],tol,[lbl,' (vertcat)']);
            S = A+B;
            assert_close(ev(S,x,allnames),VA+VB,tol,[lbl,' (plus)']);
            ncases = ncases+2;
        end
    end
end

%% Functionals (degree one)
for N = 1:2
    rng(500+N);
    [svars,dom] = spatial_space(N);
    svs = pv_cell(N);
    x = {rand_pol(svs),rand_pol(svs),rand_pol(svs)};
    specs = {   [1,1],      1
                [1,1,2],    [2,2]   };
    for k=1:size(specs,1)
        A = raw_functional(svars,dom,specs{k,1});           B = raw_functional(svars,dom,specs{k,2});
        lbl = sprintf('functional N=%d,case %d',N,k);
        VA = ev(A,x,allnames);                              VB = ev(B,x,allnames);
        [A2,B2] = common_basis(A,B);
        assert_close(ev(A2,x,allnames),VA,tol,[lbl,' (common_basis A)']);
        assert_close(ev(B2,x,allnames),VB,tol,[lbl,' (common_basis B)']);
        V = [A;B];
        assert_close(ev(V,x,allnames),[VA;VB],tol,[lbl,' (vertcat)']);
        ncases = ncases+3;
    end
end

fprintf('Test_polyopvar_common_basis_duplicates passed (%d checks).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [svars,dom] = spatial_space(N)
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
end

function svs = pv_cell(N)
pvar s1 s2 s3
svs = {s1,s2,s3};
svs = svs(1:N);
end

function p = rand_pol(svs)
p = polynomial(randn);
for i=1:numel(svs)
    p = p + randn*svs{i} + randn*svs{i}^2;
end
end

function F = raw_poly(N,Zrows,kind)
% Function with one term for each row of Zrows, as one object with Z=[Zrows] and C=[C1,C2,...], not combined: rows
% that are equal give duplicate monomials. States a,b,c (scalar) depending on s1,...,sN.
[svars,dom] = spatial_space(N);
vars = struct('in',{svars},'out',{svars});
F = [];
for k=1:numel(Zrows)
    Zrow = Zrows{k};
    d = numel(Zrow);
    if strcmp(kind,'pw'),   dims = [2*ones(d,1),ones(d,1)];     ot = [1,2];     mout = 2;
    else,                   dims = [(1+(1:d)'),ones(d,1)];      ot = [2,2];     mout = prod(dims(:,1));     end
    dp = struct('in',true(d,N),'out',true(d,N));
    Kop = rand_tensopvar(dims,vars,dom,dp,1,ot,2);
    if ~Kop.chkval,     Kop = initialize(Kop);      end
    Ck = tensopmat(Kop);
    if k==1
        F = polyopvar();
        F.C = Ck;       F.Z = Zrow;
        F.dims.out = mout;      F.depmat.out = Kop.deps.out;
    else
        n = max(size(F.Z,2),d);
        F.Z = [[F.Z,zeros(size(F.Z,1),n-size(F.Z,2))];[Zrow,zeros(1,n-d)]];
        F.C = [F.C,Ck];
    end
end
F.varname = {'a';'b';'c'};
F.dims.in = ones(3,1);
F.vars.in = svars;          F.dom.in = dom.in;          F.depmat.in = true(3,N);
F.vars.out = svars;         F.dom.out = dom.out;
F = initialize(F);
end

function F = raw_functional(svars,dom,Zcol)
% Functional with one term (random intop) for each entry of Zcol (degree one), not combined
N = numel(svars);
F = [];
for k=1:numel(Zcol)
    Aop = rand_intop([1,1],svars,dom.in,ones(1,N),2,0);
    Ck = tensopmat(Aop);
    if k==1,    F = Ck;     else,   F = [F,Ck];     end %#ok<AGROW>
end
V = polyopvar();
V.Z = Zcol(:);              V.C = F;                    V.varname = {'a';'b';'c'};
V.dims.in = ones(3,1);      V.dims.out = 1;
V.vars.in = svars;          V.vars.out = cell(1,0);
V.dom.in = dom.in;          V.dom.out = zeros(0,2);
V.depmat.in = true(3,numel(svars));     V.depmat.out = zeros(1,0);
F = initialize(V);
end

function Y = ev(P,x,allnames)
% Value of P at the states x{i} (names allnames); states of P that are not in allnames are not supported
xv = cell(numel(P.varname),1);
for k=1:numel(xv)
    xv{k} = x{strcmp(allnames,P.varname{k})};
end
Y = eval_polyopvar(P,xv);
end

function Y = sum_terms(P,x,allnames)
% Sum of the values of the individual terms (monomial j with its operators) of P
Y = 0;
for j=1:size(P.Z,1)
    T = P;
    T.Z = P.Z(j,:);
    T.C.ops = P.C.ops(:,j);
    T.C.dims.in = P.C.dims.in(j);
    T.C.depmat.in = P.C.depmat.in(j,:);
    Y = Y + ev(T,x,allnames);
end
end

function assert_close(Ynew,Yref,tol,label)
Yref = polynomial(Yref);        Ynew = polynomial(Ynew);
Yc = Yref.coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("Test_polyopvar_common_basis_duplicates (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol*max(1,max(abs(full(Yc(:)))))
    error("Test_polyopvar_common_basis_duplicates failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end
