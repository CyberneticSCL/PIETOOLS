% Test_polyopvar_dmonomials.m -- verify sosprogramming/dmonomials.m (used by several
% Examples/ files and sosprogramming/piesos_ineq.m to build monomial vectors), in two
% parts:
% (1) vectors of scalar states, for varying numbers of states and degrees. Oracle: the
%     number of monomials (nchoosek sum) and the values of the monomials, compared as
%     sorted sets (the order of the monomials in Z is not asserted). The monomials are
%     tensor products x_k1 o ... o x_kd: each factor is evaluated at its own point, so
%     that the pointwise product x(s)^d fails the test;
% (2) variables that are rows of blocks of x.C.ops: x=Top*y, linear combinations of
%     states, nonlinear rows, vector-valued states and blocks of different size, for 1
%     and 2 spatial variables, with cumulative degrees (compared as a set of values) and
%     degrees per variable (compared exactly). The monomial of degree d in a block of
%     size m is the Kronecker product of d copies of the block, of size m^d. Oracle: the
%     blocks evaluated by apply_sopvar and polynomial arithmetic, independent of the
%     'polyopvar' class; Z must have no unused monomials (empty columns of Z.C.ops), as
%     dmonomials takes variable j as vartab{j,:}; degrees per scalar row must error.

% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/22/2026: Monomials as tensor products: each factor is evaluated at its
%   own point;
% DJ (CC), 10/23/2026: Variables are the rows of blocks (part 2, merged from
%   Test_polyopvar_dmonomials_poly);

% NOTES:
% - Only correctness is tested. dmonomials.m builds each monomial by repeated '*' and
%   carries the full monomial basis seen so far (part 1).

check_scalar_states();
check_block_variables();
disp('Test_polyopvar_dmonomials passed.')


function check_scalar_states()
tol = 1e-8;

for nvars = 1:3
    for trial = 1:8
        rng(7000*nvars+trial);

        pvar s1 s2
        xs = cell(nvars,1);
        for i=1:nvars
            xs{i} = polyopvar(['x',num2str(i)],[s1,s2],[0,1]);
        end
        vartab = vertcat(xs{:});

        xvals = cell(nvars,1);
        for i=1:nvars
            xvals{i} = rand_poly([1,1],{'s1','s2'},2,4);
        end

        maxdeg = randi([0,3]);
        degs = (0:maxdeg)';
        Z = dmonomials(vartab,degs);

        % expected count = sum_{d=0}^{maxdeg} nchoosek(nvars+d-1,d)
        expected_count = 0;
        for d=0:maxdeg
            expected_count = expected_count + nchoosek(nvars+d-1,d);
        end
        if size(Z,1)~=expected_count
            error("Test_polyopvar_dmonomials failed: monomial count mismatch (nvars=%d,maxdeg=%d,trial=%d).",nvars,maxdeg,trial)
        end

        Y = apply_poly(Z,xvals);

        % values at random points, compared as sorted sets. The monomials are tensor
        % products: factor c is a function of its own copy of the variables, named
        % s1_c, s2_c (s1, s2 for degree 1), so each factor is evaluated at its own point
        pts = rand(max(maxdeg,1),2);
        xnum = zeros(nvars,size(pts,1));
        for i=1:nvars
            for c=1:size(pts,1)
                xnum(i,c) = double(subs(xvals{i},[s1;s2],pts(c,:)'));
            end
        end
        Ynum = zeros(expected_count,1);
        for k=1:expected_count
            vn = Y(k).varname(:);
            if isempty(vn)
                Ynum(k) = double(Y(k));
            else
                Ynum(k) = double(subs(Y(k),vn,cellfun(@(nm) point_of(nm,pts),vn)));
            end
        end
        % expected values: all compositions of d into nvars nonneg integer parts; the
        % factors of a monomial are ordered by state, as in dmonomials
        expected_num = [];
        for d=0:maxdeg
            combos = compositions_nonneg(d,nvars);
            for c=1:size(combos,1)
                st = repelem(1:nvars,combos(c,:));
                e = 1;
                for f=1:numel(st)
                    e = e*xnum(st(f),f);
                end
                expected_num = [expected_num; e];
            end
        end
        if numel(expected_num)~=numel(Ynum)
            error("Test_polyopvar_dmonomials failed: expected/actual count mismatch (nvars=%d,maxdeg=%d,trial=%d).",nvars,maxdeg,trial)
        end
        d1 = sort(expected_num);
        d2 = sort(Ynum);
        if max(abs(d1-d2)) > tol*max(1,max(abs(d1)))
            error("Test_polyopvar_dmonomials failed: value mismatch (nvars=%d,maxdeg=%d,trial=%d).",nvars,maxdeg,trial)
        end
    end
end
end

function check_block_variables()
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
    for icase = 1:17
        rng(100*N+icase);
        v1 = rand_pol(svs);     v2 = rand_pol(svs);
        v3 = [rand_pol(svs);rand_pol(svs);rand_pol(svs)];       % values of y3, size 3
        v4 = [rand_pol(svs);rand_pol(svs)];                     % values of y4, size 2
        deps = [];                          % variables on which each block of x depends
        switch icase
            case {1,2}
                m = icase+1;
                Top = rand_sopvar([m,1],vars,dom,1);
                x = Top*y1;                         xv = apply_sopvar(Top,v1);      bs = m;
                label = sprintf('N=%d, x=Top*y, Top %d x 1',N,m);
            case 3
                T = arrayfun(@(k) rand_sopvar([1,1],vars,dom,1),1:4,'UniformOutput',false);
                x = [T{1}*y1+T{2}*y2; T{3}*y1+T{4}*y2];
                xv = [apply_sopvar(T{1},v1)+apply_sopvar(T{2},v2); apply_sopvar(T{3},v1)+apply_sopvar(T{4},v2)];
                bs = [1,1];
                label = sprintf('N=%d, linear combination of 2 states',N);
            case 4
                T1 = rand_sopvar([1,1],vars,dom,1);
                x = [T1*y1; y1.*y2; y2];
                xv = [apply_sopvar(T1,v1); v1*v2; v2];
                bs = [1,1,1];
                label = sprintf('N=%d, nonlinear rows',N);
            case 5
                if N<2,     continue,   end
                varsS = struct('in',{{'s1','s2'}},'out',{{'s1'}});
                domS = struct('in',[0,1;0,1],'out',[0,1]);
                S = rand_sopvar([1,1],varsS,domS,1);
                x = [S*y1; y1];                     xv = [apply_sopvar(S,v1); v1];
                bs = [1,1];                         deps = logical([1 0;1 1]);
                label = 'N=2, rows of different dependence';
            case 6
                T1 = rand_sopvar([1,1],vars,dom,1);
                x = [T1*y1; 0*y1];                  xv = [apply_sopvar(T1,v1); polynomial(0)];
                bs = [1,1];
                label = sprintf('N=%d, zero row',N);
            case 7
                T1 = rand_sopvar([1,1],vars,dom,1);
                x = [T1*y1; T1*y1];                 xv = [apply_sopvar(T1,v1); apply_sopvar(T1,v1)];
                bs = [1,1];
                label = sprintf('N=%d, repeated row',N);
            case 8
                T1 = rand_sopvar([1,1],vars,dom,1);
                x = T1*y1;                          xv = apply_sopvar(T1,v1);
                bs = 1;
                label = sprintf('N=%d, single row',N);
            case 9
                x = [y1;y2];                        xv = [v1;v2];
                bs = [1,1];
                label = sprintf('N=%d, vector of states',N);
            case 10
                % A single vector-valued state of size 3 (rows of a block of C, as for u=T*x with
                % a vector state): the block is the variable of the monomials
                x = polyopvar({'y3'},3,svars,repmat([0,1],N,1),true(1,N));
                xv = v3;                            bs = 3;
                label = sprintf('N=%d, vector-valued state',N);
            case 11
                % Vector-valued state of size 3 and a scalar state: blocks of size 3 and 1
                y3 = polyopvar({'y3'},3,svars,repmat([0,1],N,1),true(1,N));
                x = [y3;y1];                        xv = [v3;v1];
                bs = [3,1];
                label = sprintf('N=%d, [y3;y1], blocks of size 3 and 1',N);
            case 12
                Top = rand_sopvar([2,1],vars,dom,1);
                x = [Top*y1;y2];                    xv = [apply_sopvar(Top,v1);v2];
                bs = [2,1];
                label = sprintf('N=%d, [Top*y1;y2], blocks of size 2 and 1',N);
            case 13
                % Scalar state, and a matrix operator on a vector-valued state
                y4 = polyopvar({'y4'},2,svars,repmat([0,1],N,1),true(1,N));
                Tm = rand_sopvar([2,2],vars,dom,1);
                x = [y2;Tm*y4];                     xv = [v2;apply_sopvar(Tm,v4)];
                bs = [1,2];
                label = sprintf('N=%d, [y2;Tm*y4], blocks of size 1 and 2',N);
            case 14
                x = polyopvar({'y4'},2,svars,repmat([0,1],N,1),true(1,N));
                xv = v4;                            bs = 2;
                label = sprintf('N=%d, vector-valued state of size 2',N);
            case 15
                % A block of size 2 and a nonlinear (pointwise product) row
                Top = rand_sopvar([2,1],vars,dom,1);
                x = [Top*y1;y1.*y2];                xv = [apply_sopvar(Top,v1);v1*v2];
                bs = [2,1];
                label = sprintf('N=%d, [Top*y1;y1.*y2], block of size 2 and nonlinear row',N);
            case 16
                if N<2,     continue,   end
                varsS = struct('in',{{'s1','s2'}},'out',{{'s1'}});
                domS = struct('in',[0,1;0,1],'out',[0,1]);
                S = rand_sopvar([2,1],varsS,domS,1);
                x = [S*y1; y1];                     xv = [apply_sopvar(S,v1); v1];
                bs = [2,1];                         deps = logical([1 0;1 1]);
                label = 'N=2, block of size 2 and row of different dependence';
            case 17
                % Two vector-valued states, of size 3 and 2
                y3 = polyopvar({'y3'},3,svars,repmat([0,1],N,1),true(1,N));
                y4 = polyopvar({'y4'},2,svars,repmat([0,1],N,1),true(1,N));
                x = [y3;y4];                        xv = [v3;v4];
                bs = [3,2];
                label = sprintf('N=%d, [y3;y4], blocks of size 3 and 2',N);
        end
        x0 = x;     if ~x0.chkval,  x0 = initialize(x0);    end
        nr = size(x0,1);                    % number of rows
        nv = numel(bs);                     % number of blocks = variables of the monomials
        if isempty(deps),   deps = true(nv,N);      end
        if nr~=numel(xv) || nr~=sum(bs) || ~isequal(reshape(x0.dims.out,1,[]),bs)
            error("Test_polyopvar_dmonomials (%s): size(x,1)=%d, dims.out=%s, expected %d, %s.",label,nr,mat2str(x0.dims.out),numel(xv),mat2str(bs))
        end
        % values of the states by name, with the vector-valued states as columns
        vals = struct('y1',v1,'y2',v2,'y3',v3,'y4',v4);
        xb = split_rows(xv,bs);             % values of the blocks, the factors of the monomials
        for dmax = 1:3
            if N==2 && sum(bs)*dmax>6 && any(icase==[4,15]),   continue,   end     % skip: nonlinear rows, high degree
            % (a) cumulative degrees, as a set of values
            degs = (0:dmax)';
            Z = dmonomials(x,degs);
            nZ_exp = 0;     sizes_exp = [];
            for d=0:dmax
                combos = compositions_nonneg(d,nv);
                nZ_exp = nZ_exp + size(combos,1);
                sizes_exp = [sizes_exp;prod(repmat(bs(:)',size(combos,1),1).^combos,2)];       %#ok<AGROW>
            end
            if ~Z.chkval,   Z = initialize(Z);      end
            assert_no_unused(Z,sprintf('%s, d=%d',label,dmax));
            if numel(Z.dims.out)~=nZ_exp || size(Z,1)~=sum(sizes_exp) || ~isequal(sort(Z.dims.out(:)),sort(sizes_exp))
                error("Test_polyopvar_dmonomials failed (%s, d=%d): %d monomials of %d rows, expected %d of %d.",label,dmax,numel(Z.dims.out),size(Z,1),nZ_exp,sum(sizes_exp))
            end
            Yz = eval_Z(Z,vals);
            Ye = [];
            for d=0:dmax
                combos = compositions_nonneg(d,nv);
                Ye = [Ye; expected_tensor(xb,combos,deps,svars)];       %#ok<AGROW>
            end
            % random values for all variables, including the copies s1_c of the factors
            names = unique([Yz.varname(:);Ye.varname(:)]);
            mp = containers.Map('KeyType','char','ValueType','double');
            for k=1:numel(names),   mp(names{k}) = rand;    end
            a = sort(values_at(Yz,mp));
            b = sort(values_at(Ye,mp));
            assert_close(a,b,tol,sprintf('%s, cumulative degrees 0:%d',label,dmax));
            ncases = ncases+1;
        end
        % (b) degrees of each variable, compared exactly, in the given order
        D = [2,zeros(1,nv-1); ones(1,nv); zeros(1,nv-1),1; 1,zeros(1,nv-1)];
        if nv==1,   D = [2;1;0;3];  end
        if nv==3,   D(end+1,:) = [1,0,2];   end
        Z = dmonomials(x,D);
        if ~Z.chkval,   Z = initialize(Z);      end
        assert_no_unused(Z,sprintf('%s, degree array',label));
        if numel(Z.dims.out)~=size(D,1)
            error("Test_polyopvar_dmonomials failed (%s): %d monomials, expected %d.",label,numel(Z.dims.out),size(D,1))
        end
        Yz = eval_Z(Z,vals);
        Ye = expected_tensor(xb,D,deps,svars);
        assert_poly_equal(Yz,Ye,tol,sprintf('%s, degree array',label));
        ncases = ncases+1;
        % (c) degrees per scalar row of a vector-valued variable are not degrees per variable
        if nr>nv
            Dr = zeros(1,nr);   Dr(1) = 1;
            assert_error(@() dmonomials(x,[Dr;Dr]),'Degrees must be specified',sprintf('%s, degrees per scalar row',label));
            ncases = ncases+1;
        end
    end
end
end

function v = point_of(nm,pts)
% Value of the variable s1 or s2, or its copy c, s1_c (as named by apply_op), at the
% independent points pts(c,:)
tok = regexp(nm,'^s(\d)_?(\d*)$','tokens','once');
c = 1;
if ~isempty(tok{2}),    c = str2double(tok{2});     end
v = pts(c,str2double(tok{1}));
end


function combos = compositions_nonneg(d,n)
% All n-tuples of nonnegative integers summing to d (rows of an array).
if n==1
    combos = d;
    return
end
combos = zeros(0,n);
for k=0:d
    sub = compositions_nonneg(d-k,n-1);
    combos = [combos; [k*ones(size(sub,1),1), sub]];
end
end

function p = rand_pol(svs)
p = polynomial(randn);
for i=1:numel(svs)
    p = p + randn*svs{i} + randn*svs{i}^2;
end
end


function xb = split_rows(xv,bs)
% Cell with the blocks of rows of the polynomial column xv, of sizes bs
xb = cell(1,numel(bs));
r0 = 0;
for j=1:numel(bs)
    xb{j} = xv(r0+(1:bs(j)));
    r0 = r0+bs(j);
end
end


function Y = eval_Z(Z,vals)
% Z applied to the states (by name), as a column of polynomials
xc = cellfun(@(nm) vals.(nm),reshape(Z.varname,[],1),'UniformOutput',false);
Y = polynomial(apply_poly(Z,xc));
end


function c = kron_poly(a,b)
% Kronecker product of the polynomial columns a and b: c((i-1)*numel(b)+j) = a(i)*b(j)
nb = numel(b);
c = polynomial(zeros(numel(a)*nb,1));
for i=1:numel(a)
    for j=1:nb
        c((i-1)*nb+j) = a(i)*b(j);
    end
end
end


function E = expected_tensor(xb,degmat,deps,svars)
% Monomial i, the Kronecker product of its factors xb{j}, j repeated degmat(i,j) times
% in the order of the variables, stacked as a column of polynomials. As apply_op names
% them for a monomial of degree >1, the variable s_v that factor f depends on (deps(j,v))
% is renamed to s_v_c, c the number of factors up to f depending on s_v
nZ = size(degmat,1);
E = [];
for i=1:nZ
    rows = repelem(1:size(degmat,2),degmat(i,:));
    cnt = zeros(1,numel(svars));
    e = polynomial(1);
    for f=1:numel(rows)
        xf = xb{rows(f)};
        cnt = cnt+deps(rows(f),:);
        if numel(rows)>1
            for v=find(deps(rows(f),:))
                xf.varname(strcmp(xf.varname,svars{v})) = {[svars{v},'_',num2str(cnt(v))]};
            end
        end
        e = kron_poly(e,xf);
    end
    E = [E;e];          %#ok<AGROW>
end
end


function v = values_at(P,mp)
% Values of the polynomials in the column P at the values mp(name) of the variables
v = zeros(numel(P),1);
for k=1:numel(P)
    vn = P(k).varname(:);
    if isempty(vn)
        v(k) = double(P(k));
    else
        vals = cell2mat(values(mp,reshape(vn,1,[])));
        v(k) = double(subs(P(k),vn,vals(:)));
    end
end
end


function assert_no_unused(Z,label)
% Every monomial (column of Z.C.ops) must be acted on by some block: the variables are
% extracted as the nonempty terms of their rows of blocks, vartab{j,:}, not as vartab(rows,:),
% which keeps all monomials of vartab and so leaves unused monomials in the products
unused = find(~any(~cellfun(@isempty,Z.C.ops),1));
if ~isempty(unused)
    error("Test_polyopvar_dmonomials failed (%s): monomial(s) %s of Z are not acted on by any block.",label,mat2str(unused))
end
end


function assert_close(a,b,tol,label)
if numel(a)~=numel(b) || max(abs(a-b))>tol*max(1,max(abs(b)))
    error("Test_polyopvar_dmonomials failed (%s): values differ (max error %.3e).",label,max(abs(a-b)))
end
end


function assert_error(f,needle,label)
try
    f();
catch err
    if ~contains(err.message,needle)
        error('Test_polyopvar_dmonomials failed (%s): wrong error "%s".',label,err.message)
    end
    return
end
error('Test_polyopvar_dmonomials failed (%s): expected an error containing "%s".',label,needle)
end


function assert_poly_equal(Y,E,tol,label)
if numel(Y)~=numel(E)
    error("Test_polyopvar_dmonomials failed (%s): %d entries, expected %d.",label,numel(Y),numel(E))
end
dd = polynomial(Y-E);
err = max([0;abs(full(dd.coeff(:)))]);
if err>tol
    error("Test_polyopvar_dmonomials failed (%s): residual %.3e.",label,err)
end
% reference must not be zero everywhere (vacuous test)
Ec = polynomial(E).coeff;
if isempty(Ec) || max(abs(full(Ec(:))))<1e-6
    error("Test_polyopvar_dmonomials (%s): reference is zero; test would be vacuous.",label)
end
end
