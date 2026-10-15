% Test_polyopvar_dmonomials_poly.m -- dmonomials(x,degs) for x whose rows
% are polynomials in the states, e.g. x=Top*y for a 'sopvar' Top, rather
% than a vector of states:
%   (i)   x=Top*y, Top m x 1 (m=2,3): a single block of m rows in C (all
%         rows must be used);
%   (ii)  rows that are linear combinations of two states,
%         [T1*y1+T2*y2; T3*y1+T4*y2];
%   (iii) nonlinear rows, [T1*y1; y1.*y2; y2];
%   (iv)  rows that are functions of different variables, [S*y1; y1], S
%         integrating s2;
%   (v)   a zero row, a repeated row, and a single row;
%   (vi)  a vector of states (the columns of degs must follow the rows
%         of x);
% for 1 and 2 spatial variables, with the degrees specified as cumulative
% degrees (0:d)' (compared as a set of values at random points, since the
% order is not specified) and as an n x p array of degrees of the p rows
% (compared exactly, including the order).
%
% Verification is semantic: the rows x_j are evaluated by applying the
% sopvar operators to polynomial states (apply_sopvar, independent of the
% 'polyopvar' class), and the monomials Z_i = prod_j x_j^degmat(i,j)
% formed by polynomial arithmetic, then compared to Z applied to the same
% states (apply_poly).

% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;
pvar s1 s2
svall = {s1,s2};
ncases = 0;

for N = 1:2
    svars = {'s1','s2'};        svars = svars(1:N);
    svs = svall(1:N);
    vars = struct('in',{svars},'out',{svars});
    dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
    pvs = svs{1};
    for i=2:N,  pvs = [pvs,svs{i}];     end
    y1 = polyopvar('y1',pvs,[0,1]);
    y2 = polyopvar('y2',pvs,[0,1]);
    for icase = 1:9
        rng(100*N+icase);
        v1 = rand_pol(svs);     v2 = rand_pol(svs);
        switch icase
            case {1,2}
                m = icase+1;
                Top = rand_sopvar([m,1],vars,dom,1);
                x = Top*y1;                         xv = apply_sopvar(Top,v1);
                label = sprintf('N=%d, x=Top*y, Top %d x 1',N,m);
            case 3
                T = arrayfun(@(k) rand_sopvar([1,1],vars,dom,1),1:4,'UniformOutput',false);
                x = [T{1}*y1+T{2}*y2; T{3}*y1+T{4}*y2];
                xv = [apply_sopvar(T{1},v1)+apply_sopvar(T{2},v2); apply_sopvar(T{3},v1)+apply_sopvar(T{4},v2)];
                label = sprintf('N=%d, linear combination of 2 states',N);
            case 4
                T1 = rand_sopvar([1,1],vars,dom,1);
                x = [T1*y1; y1.*y2; y2];
                xv = [apply_sopvar(T1,v1); v1*v2; v2];
                label = sprintf('N=%d, nonlinear rows',N);
            case 5
                if N<2,     continue,   end
                varsS = struct('in',{{'s1','s2'}},'out',{{'s1'}});
                domS = struct('in',[0,1;0,1],'out',[0,1]);
                S = rand_sopvar([1,1],varsS,domS,1);
                x = [S*y1; y1];                     xv = [apply_sopvar(S,v1); v1];
                label = 'N=2, rows of different dependence';
            case 6
                T1 = rand_sopvar([1,1],vars,dom,1);
                x = [T1*y1; 0*y1];                  xv = [apply_sopvar(T1,v1); polynomial(0)];
                label = sprintf('N=%d, zero row',N);
            case 7
                T1 = rand_sopvar([1,1],vars,dom,1);
                x = [T1*y1; T1*y1];                 xv = [apply_sopvar(T1,v1); apply_sopvar(T1,v1)];
                label = sprintf('N=%d, repeated row',N);
            case 8
                T1 = rand_sopvar([1,1],vars,dom,1);
                x = T1*y1;                          xv = apply_sopvar(T1,v1);
                label = sprintf('N=%d, single row',N);
            case 9
                x = [y1;y2];                        xv = [v1;v2];
                label = sprintf('N=%d, vector of states',N);
        end
        nr = size(x,1);
        if nr~=numel(xv)
            error("Test_polyopvar_dmonomials_poly (%s): size(x,1)=%d, expected %d.",label,nr,numel(xv))
        end
        vals = struct('y1',v1,'y2',v2);
        for dmax = 1:3
            if N==2 && nr*dmax>6 && icase==4,   continue,   end     % skip degree 2*dmax in 3 states
            % (a) cumulative degrees, as a set of values
            degs = (0:dmax)';
            Z = dmonomials(x,degs);
            nZ_exp = 0;
            for d=0:dmax,   nZ_exp = nZ_exp + nchoosek(nr+d-1,d);   end
            if size(Z,1)~=nZ_exp
                error("Test_polyopvar_dmonomials_poly failed (%s, d=%d): %d monomials, expected %d.",label,dmax,size(Z,1),nZ_exp)
            end
            Yz = eval_Z(Z,vals);
            Ye = [];
            for d=0:dmax
                combos = compositions_nonneg(d,nr);
                Ye = [Ye; expected_monomials(xv,combos)];       %#ok<AGROW>
            end
            pt = rand(N,1);
            a = sort(double(subs(Yz,svs_col(svs),pt)));
            b = sort(double(subs(Ye,svs_col(svs),pt)));
            assert_close(a,b,tol,sprintf('%s, cumulative degrees 0:%d',label,dmax));
            ncases = ncases+1;
        end
        % (b) degrees of each row, compared exactly, in the given order
        D = [2,zeros(1,nr-1); ones(1,nr); zeros(1,nr-1),1; 1,zeros(1,nr-1)];
        if nr==1,   D = [2;1;0;3];  end
        if nr==3,   D(end+1,:) = [1,0,2];   end
        Z = dmonomials(x,D);
        if size(Z,1)~=size(D,1)
            error("Test_polyopvar_dmonomials_poly failed (%s): %d monomials, expected %d.",label,size(Z,1),size(D,1))
        end
        Yz = eval_Z(Z,vals);
        Ye = expected_monomials(xv,D);
        assert_poly_equal(Yz,Ye,tol,sprintf('%s, degree array',label));
        ncases = ncases+1;
    end
end

fprintf('Test_polyopvar_dmonomials_poly passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function p = rand_pol(svs)
p = polynomial(randn);
for i=1:numel(svs)
    p = p + randn*svs{i} + randn*svs{i}^2;
end
end

function pv = svs_col(svs)
pv = svs{1};
for i=2:numel(svs),     pv = [pv;svs{i}];       end
end

function Y = eval_Z(Z,vals)
% Z applied to the states (by name), as a column of polynomials
xc = cellfun(@(nm) vals.(nm),reshape(Z.varname,[],1),'UniformOutput',false);
Y = polynomial(apply_poly(Z,xc));
end

function E = expected_monomials(xv,degmat)
% Rows prod_j xv(j)^degmat(i,j) as a column of polynomials
nZ = size(degmat,1);
E = polynomial(zeros(nZ,1));
for i=1:nZ
    e = polynomial(1);
    for j=1:size(degmat,2)
        if degmat(i,j)>0,   e = e*xv(j)^degmat(i,j);    end
    end
    E(i) = e;
end
end

function combos = compositions_nonneg(d,n)
% All n-tuples of nonnegative integers summing to d (rows of an array)
if n==1
    combos = d;
    return
end
combos = zeros(0,n);
for k=0:d
    sub = compositions_nonneg(d-k,n-1);
    combos = [combos; [k*ones(size(sub,1),1), sub]];        %#ok<AGROW>
end
end

function assert_close(a,b,tol,label)
if numel(a)~=numel(b) || max(abs(a-b))>tol*max(1,max(abs(b)))
    error("Test_polyopvar_dmonomials_poly failed (%s): values differ (max error %.3e).",label,max(abs(a-b)))
end
end

function assert_poly_equal(Y,E,tol,label)
dd = polynomial(Y-E);
err = max([0;abs(full(dd.coeff(:)))]);
if err>tol
    error("Test_polyopvar_dmonomials_poly failed (%s): residual %.3e.",label,err)
end
% reference must not be zero everywhere (vacuous test)
Ec = polynomial(E).coeff;
if isempty(Ec) || max(abs(full(Ec(:))))<1e-6
    error("Test_polyopvar_dmonomials_poly (%s): reference is zero; test would be vacuous.",label)
end
end
