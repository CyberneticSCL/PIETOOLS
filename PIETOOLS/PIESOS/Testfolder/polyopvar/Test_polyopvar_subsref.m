% Test_polyopvar_subsref.m -- verify @polyopvar/subsref.m, in two parts. A 'polyopvar' p is a
% vector-valued (matdim(2)=1) function p(x) = C*Z(x), whose blocks are the rows of blocks of C.ops:
% (1) '()': p(i,j) is element (i,j) of the matrix-valued polynomial, so p(i), p(i,1), p(i,:) select
%     the scalar row i, and p(:), p([i1,i2]) several rows; the column index can only be 1 or ':';
% (2) '{}': p{i,j} is block (i,j) of the matrix-valued polynomial, so p{i}, p{i,1}, p{i,:} are the
%     polynomial with the coefficients p.C{i,:}, acting on the nonempty terms (monomials) of block row i
%     (the monomials of Z without nonempty block in row i are dropped); p{:,1}, p{[i1,i2],1} several
%     blocks; the column index can only be 1 or ':'. Checked for p with diagonal blocks, [y1;y2],
%     [T1*y1;y1.*y2;y2], vector-valued blocks, y3, [y3;y1], [Top*y1;y2], rows with several terms, and
%     1 and 2 spatial variables. The values of the rows and blocks, xv, xb{j}, come from the operators
%     and polynomials directly (apply_sopvar), independent of apply_poly(p) and of the class under test.

% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: '{}' with ':' and numeric column vectors, diagonal and
%   vector-valued blocks, multiple rows, linear index over blocks; f{:,k} for rows
%   without a block; '()' part merged from Test_polyopvar_cat_subsref (vertcat is in
%   Test_polyopvar_vertcat);
% DJ (CC), 11/08/2026: New meaning of the indices, p(i,j) is element (i,j) and p{i,j} block (i,j) of
%   the matrix-valued polynomial (only vector-valued, so j=1): p(i), p(i,1), p{i}, p{i,1}; p{i,k}
%   with k a monomial and p{:,k} are no longer supported (were: terms of the blocks C.ops{i,k});

% NOTES:
% - Only vector-valued polynomials are supported: a column index other than 1 or ':' is an error,
%   for '()' and '{}'.
% - p{i,j} must return a 'tensopmat' coefficient (not the bare tensopvar/intop that C.ops{i,j}
%   unwraps to), else out.dims.out / out.depmat.out fail.
% - 'end' evaluates to 1 for a 'polyopvar' (MATLAB's size of the scalar object): p{end} is p{1}, p(end)
%   is p(1). Not supported, and not tested except for this value.

check_paren_subsref();
check_curly_subsref();
disp('Test_polyopvar_subsref passed.')


function check_paren_subsref()
tol = 1e-8;

pvar s1 s2
x = polyopvar('x',[s1,s2],[0,1]);
y = polyopvar('y',[s1,s2],[0,1]);

xv = {rand_poly([1,1],{'s1','s2'},2,5)};
yv = {rand_poly([1,1],{'s1','s2'},2,5)};

v = [x;y];

% '()' subsref, rows, with the column index 1 or ':', or a single index
for vk = {v(2,:),v(2,1),v(2),v(2,':')}
    Vv2 = apply_poly(vk{1},{xv{1},yv{1}});
    errv2 = cleanpoly(Vv2-yv{1},1e-10);
    if max(max(abs(errv2.C)))>tol*max(1,max(max(abs(yv{1}.C))))
        error("Test_polyopvar_subsref failed: row-subsref value mismatch.")
    end
end

% a scalar-valued function has the single element f(1,1), also if it has several terms
f = x + x.*x;
Vf = apply_poly(f(1,1),xv);
errf = cleanpoly(Vf-(xv{1}+xv{1}*xv{1}),1e-10);
if max(max(abs(errf.C)))>tol*max(1,max(max(abs(xv{1}.C))))
    error("Test_polyopvar_subsref failed: f(1,1) is not f.")
end
if ~same_terms(f(1,1),f) || ~same_terms(f(1),f) || ~same_terms(f(:),f)
    error("Test_polyopvar_subsref failed: f(1,1), f(1) and f(:) must retain all terms.")
end

% a column index other than 1 (only vector-valued polynomials are supported)
assert_error(@() f(1,2),'column index must be 1',"f(1,2)");
assert_error(@() v(2,2),'column index must be 1',"v(2,2)");
assert_error(@() v(1,[1,2]),'column index must be 1',"v(1,[1,2])");
end

function check_curly_subsref()
tol = 1e-8;

pvar s1 s2
x = polyopvar('x',[s1,s2],[0,1]);
y = polyopvar('y',[s1,s2],[0,1]);
f = x + x.*x;       % 2 distinct monomial terms: bare x, and x^2
v = [f;y];          % block row 1: f's 2 terms; block row 2: y's single (different) term

[nr,nc] = size(v.C.ops);
if nr~=2 || nc<3
    error("Test_polyopvar_subsref failed: test setup expected 2 rows of blocks and >=3 monomial columns.")
end

for trial=1:5

xv = rand_poly([1,1],{'s1','s2'},2,5);
yv = rand_poly([1,1],{'s1','s2'},2,5);
Yfull = apply_poly(v,{xv,yv});

% block i is the ith entry: the polynomial with the coefficients v.C{i,:} on the nonempty terms
for i=1:nr
    vi = v{i,1};
    if ~isa(vi,'polyopvar')
        error("Trial %d: v{%d,1} did not return a 'polyopvar' object.",trial,i)
    end
    pj = find(~cellfun(@isempty,v.C.ops(i,:)));
    if ~isequal(vi.Z,v.Z(pj,:))
        error("Trial %d: v{%d,1}.Z does not match the monomials v.Z of the nonempty terms.",trial,i)
    end
    erri = cleanpoly(apply_poly(vi,{xv,yv})-Yfull(i),1e-10);
    if max(max(abs(erri.C)))>tol*max(1,max(max(abs(Yfull(i).C))))
        error("Trial %d: v{%d,1} does not represent row %d of v.",trial,i,i)
    end
end

end

%%%% Linear indexing must match the corresponding (row,1) pair
for k=1:nr
    vk = v{k};      vk1 = v{k,1};       vkc = v{k,:};
    if ~isequal(vk.Z,vk1.Z) || ~isequal(vk.C.ops,vk1.C.ops) || ~isequal(vkc.Z,vk1.Z) || ~isequal(vkc.C.ops,vk1.C.ops)
        error("Test_polyopvar_subsref failed: v{%d}, v{%d,1} and v{%d,:} differ.",k,k,k)
    end
end

%%%% Out-of-range indices, and column indices other than 1, must error
assert_error(@() v{nr+1},'Row indices cannot exceed',"row beyond the blocks, v{nr+1}");
assert_error(@() v{nr+1,1},'Row indices cannot exceed',"row beyond the blocks, v{nr+1,1}");
assert_error(@() v{1,2},'column index must be 1',"v{1,2}");
assert_error(@() v{1,nc},'column index must be 1',"v{1,nc}");
assert_error(@() v{1,0},'column index must be 1',"v{1,0}");
assert_error(@() v{2,[1,1]},'column index must be 1',"v{2,[1,1]}");

%%%% ':' and numeric row vectors, diagonal and vector-valued blocks
ncases = 0;
svall = {s1,s2};
for N=1:2
    svars = {'s1','s2'};        svars = svars(1:N);
    svs = svall(1:N);
    vars = struct('in',{svars},'out',{svars});
    dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
    pvs = svs{1};
    for k=2:N,  pvs = [pvs,svs{k}];     end
    ya = polyopvar('y1',pvs,[0,1]);
    yb = polyopvar('y2',pvs,[0,1]);
    y3 = polyopvar({'y3'},3,svars,repmat([0,1],N,1),true(1,N));
    y4 = polyopvar({'y4'},2,svars,repmat([0,1],N,1),true(1,N));
    for icase=1:12
        rng(100*N+icase);
        v1 = rand_pol(svs);     v2 = rand_pol(svs);
        v3 = [rand_pol(svs);rand_pol(svs);rand_pol(svs)];       % values of y3, size 3
        v4 = [rand_pol(svs);rand_pol(svs)];                     % values of y4, size 2
        Top = rand_sopvar([2,1],vars,dom,1);
        T1 = rand_sopvar([1,1],vars,dom,1);
        switch icase
            case 1
                x = [ya;yb];                    xv = [v1;v2];                       bs = [1,1];
                label = sprintf('N=%d, [y1;y2]',N);
            case 2
                x = [y3;ya];                    xv = [v3;v1];                       bs = [3,1];
                label = sprintf('N=%d, [y3;y1]',N);
            case 3
                x = [Top*ya;yb];                xv = [apply_sopvar(Top,v1);v2];     bs = [2,1];
                label = sprintf('N=%d, [Top*y1;y2]',N);
            case 4
                x = [T1*ya;ya.*yb;yb];          xv = [apply_sopvar(T1,v1);v1*v2;v2];    bs = [1,1,1];
                label = sprintf('N=%d, [T1*y1;y1.*y2;y2]',N);
            case 5
                x = y3;                         xv = v3;                            bs = 3;
                label = sprintf('N=%d, y3, a single vector-valued block',N);
            case 6
                x = Top*ya;                     xv = apply_sopvar(Top,v1);          bs = 2;
                label = sprintf('N=%d, Top*y1, a single block of 2 rows',N);
            case 7
                x = [y3;y4];                    xv = [v3;v4];                       bs = [3,2];
                label = sprintf('N=%d, [y3;y4]',N);
            case 8
                Tm = rand_sopvar([2,2],vars,dom,1);
                x = [yb;Tm*y4];                 xv = [v2;apply_sopvar(Tm,v4)];      bs = [1,2];
                label = sprintf('N=%d, [y2;Tm*y4]',N);
            case 9
                x = [T1*ya+ya.*yb;yb];          xv = [apply_sopvar(T1,v1)+v1*v2;v2];    bs = [1,1];
                label = sprintf('N=%d, [T1*y1+y1.*y2;y2], a row with 2 terms',N);
            case 10
                x = [Top*ya;ya.*yb];            xv = [apply_sopvar(Top,v1);v1*v2];  bs = [2,1];
                label = sprintf('N=%d, [Top*y1;y1.*y2]',N);
            case 11
                if N<2,     continue,   end
                varsS = struct('in',{{'s1','s2'}},'out',{{'s1'}});
                domS = struct('in',[0,1;0,1],'out',[0,1]);
                S = rand_sopvar([1,1],varsS,domS,1);
                x = [S*ya;ya];                  xv = [apply_sopvar(S,v1);v1];       bs = [1,1];
                label = 'N=2, [S*y1;y1], rows of different dependence';
            case 12
                x = [y3;ya;y4];                 xv = [v3;v1;v4];                    bs = [3,1,2];
                label = sprintf('N=%d, [y3;y1;y4]',N);
        end
        v = x;
        if ~v.chkval,   v = initialize(v);  end
        nb = numel(bs);
        [nbr,nbc] = size(v.C.ops);          % blocks of C.ops, not scalar rows/columns of C
        if nbr~=nb || ~isequal(reshape(v.dims.out,1,[]),bs) || size(v,1)~=sum(bs)
            error("Test_polyopvar_subsref (%s): unexpected block structure of the test function.",label)
        end
        pop = ~cellfun(@isempty,v.C.ops);   % populated blocks
        vals = struct('y1',v1,'y2',v2,'y3',v3,'y4',v4);
        xc = cellfun(@(nm) vals.(nm),reshape(v.varname,[],1),'UniformOutput',false);
        xb = split_rows(xv,bs);
        % apply_poly of the whole function must agree with the values of the blocks
        assert_poly_equal(apply_poly(v,xc),xv,tol,[label ', apply_poly(v)']);

        % v(k), v(k,1), v(k,:): the scalar row k of the matrix-valued polynomial v, with all monomials
        nrv = sum(bs);
        for k=1:nrv
            for vk = {v(k),v(k,1),v(k,:)}
                assert_poly_equal(apply_poly(vk{1},xc),xv(k),tol,sprintf('%s, v(%d)',label,k));
                if size(vk{1},1)~=1 || size(vk{1}.Z,1)~=size(v.Z,1)
                    error("Test_polyopvar_subsref failed (%s): v(%d) is not a single row with all monomials.",label,k)
                end
            end
            if ~same_terms(v(k),v(k,:)) || ~same_terms(v(k,1),v(k,:))
                error("Test_polyopvar_subsref failed (%s): v(%d), v(%d,1) and v(%d,:) differ.",label,k,k,k)
            end
            ncases = ncases+1;
        end
        % v(:), v(:,1), v(:,:) are v, several rows in the given order
        for vk = {v(:),v(:,1),v(:,:)}
            assert_poly_equal(apply_poly(vk{1},xc),xv,tol,[label ', v(:)']);
            if ~same_terms(vk{1},v)
                error("Test_polyopvar_subsref failed (%s): v(:) is not v.",label)
            end
        end
        if nrv>=2
            r = [nrv,1];
            assert_poly_equal(apply_poly(v(r),xc),xv(r),tol,[label ', v([nrv,1])']);
            assert_poly_equal(apply_poly(v(r,1),xc),xv(r),tol,[label ', v([nrv,1],1)']);
            r = 2:nrv;
            assert_poly_equal(apply_poly(v(r),xc),xv(r),tol,[label ', v(2:end)']);
        end
        assert_error(@() v(nrv+1),'exceed',[label ', row beyond the rows']);
        assert_error(@() v(0),'positive',[label ', row 0']);
        assert_error(@() v(1,2),'column index must be 1',[label ', column beyond the single column']);
        assert_error(@() v(1,[1,1]),'column index must be 1',[label ', repeated column']);
        assert_error(@() v(1,1,1),'At most two indices',[label ', three indices']);
        ncases = ncases+6;

        % v{j}, v{j,1}, v{j,:}: block j of v, the polynomial with the coefficients of row of blocks j,
        % on its nonempty terms and their monomials
        r0 = 0;
        for j=1:nb
            vj = v{j,1};
            lbl = sprintf('%s, v{%d,1}',label,j);
            if ~isa(vj,'polyopvar') || size(vj,1)~=bs(j)
                error("Test_polyopvar_subsref failed (%s): wrong class or size.",lbl)
            elseif ~isequal(size(vj.C.ops),[1,nnz(pop(j,:))]) || any(cellfun(@isempty,vj.C.ops(:)))
                error("Test_polyopvar_subsref failed (%s): not the nonempty terms (%s).",lbl,mat2str(size(vj.C.ops)))
            elseif ~isequal(full(vj.Z),full(v.Z(pop(j,:),:)))
                error("Test_polyopvar_subsref failed (%s): monomials do not match those of the nonempty terms.",lbl)
            end
            assert_poly_equal(apply_poly(vj,xc),xb{j},tol,lbl);
            % the linear index, and ':', must give the same as the column index 1
            if ~same_terms(vj,v{j}) || ~same_terms(vj,v{j,:})
                error("Test_polyopvar_subsref failed (%s): v{j} or v{j,:} differs from v{j,1}.",lbl)
            end
            assert_poly_equal(apply_poly(v{j},xc),xb{j},tol,sprintf('%s, v{%d}',label,j));
            % for a scalar block, v(j) is the block
            rows_j = r0+(1:bs(j));      r0 = r0+bs(j);
            assert_poly_equal(apply_poly(v(rows_j),xc),xb{j},tol,sprintf('%s, v(rows of block %d)',label,j));
            ncases = ncases+1;
        end

        % stacking v{j,1} over all j, and v{:,1}, v{:,:}, reproduce v
        Vs = v{1,1};
        for j=2:nb,     Vs = [Vs;v{j,1}];    end                %#ok<AGROW>
        assert_poly_equal(apply_poly(Vs,xc),xv,tol,[label ', [v{1};v{2};...]']);
        % (assigned first: in a cell array, MATLAB would expect one output per selected block)
        Va1 = v{:,1};   Va2 = v{:,:};   Va3 = v{1:nb,1};
        for Va = {Va1,Va2,Va3}
            assert_poly_equal(apply_poly(Va{1},xc),xv,tol,[label ', v{:,1}']);
            if ~isequal(full(Va{1}.Z),full(v.Z))
                error("Test_polyopvar_subsref failed (%s): v{:,1} does not retain all monomials.",label)
            end
        end
        % several blocks at once, in the given order
        if nb>=2
            for rows = {[1,nb],[nb,1],2:nb}
                r = rows{1};
                vr = v{r,1};
                assert_poly_equal(apply_poly(vr,xc),vertcat(xb{r}),tol,sprintf('%s, v{%s,1}',label,mat2str(r)));
                vr = v{r,:};
                assert_poly_equal(apply_poly(vr,xc),vertcat(xb{r}),tol,sprintf('%s, v{%s,:}',label,mat2str(r)));
            end
        end

        % indices beyond the blocks (the scalar dimensions size(v.C) are larger for vector-valued blocks),
        % and only the single column of blocks
        assert_error(@() v{nb+1,1},'Row indices cannot exceed',[label ', row beyond the blocks']);
        assert_error(@() v{nb+1},'Row indices cannot exceed',[label ', linear index beyond the blocks']);
        assert_error(@() v{1,2},'column index must be 1',[label ', column beyond the single column']);
        if nbc>1        % else v{1,nbc} is v{1,1}
            assert_error(@() v{1,nbc},'column index must be 1',[label ', column of monomials']);
        end
        ncases = ncases+5;
    end
end

%%%% Empty selections, and the indices that remain unsupported
vv = [ya;yb];
assert_error(@() vv{':'},'Linear indexing',"linear ':'");
assert_error(@() vv{[1,2]},'Linear indexing',"linear index vector");
assert_error(@() vv{1,1,1},'At most two indices',"three indices");
assert_error(@() vv{[],1},'no nonempty term',"empty row selection");
assert_error(@() vv{1,[]},'column index must be 1',"empty column selection");
% subsequent references are applied to the extracted object
if ~isequal(full(vv{2}.Z),2) || ~isequal(vv{2,:}.dims.out,1) || ~isequal(full(vv{2,1}.Z),2)
    error("Test_polyopvar_subsref failed: chained reference after vv{2}.")
end
% 'end' evaluates to 1, see NOTES
if ~same_terms(vv{end},vv{1}) || ~same_terms(vv{1,end},vv{1,1}) || ~same_terms(vv(end),vv(1))
    error("Test_polyopvar_subsref failed: 'end' no longer evaluates to 1, update the NOTES.")
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


function tf = same_terms(f,g)
% Same monomials and the same pattern of blocks in the 'polyopvar' objects f and g
tf = isequal(size(f),size(g)) && isequal(full(f.Z),full(g.Z)) && isequal(size(f.C.ops),size(g.C.ops)) ...
        && isequal(cellfun(@isempty,f.C.ops),cellfun(@isempty,g.C.ops));
end


function assert_poly_equal(Y,E,tol,label)
if numel(Y)~=numel(E)
    error("Test_polyopvar_subsref failed (%s): %d entries, expected %d.",label,numel(Y),numel(E))
end
dd = polynomial(Y-E);
err = max([0;abs(full(dd.coeff(:)))]);
if err>tol
    error("Test_polyopvar_subsref failed (%s): residual %.3e.",label,err)
end
% reference must not be zero everywhere (vacuous test)
Ec = polynomial(E).coeff;
if isempty(Ec) || max(abs(full(Ec(:))))<1e-6
    error("Test_polyopvar_subsref (%s): reference is zero; test would be vacuous.",label)
end
end


function assert_error(f,needle,label)
try
    f();
catch err
    if ~contains(err.message,needle)
        error('Test_polyopvar_subsref failed (%s): wrong error "%s".',label,err.message)
    end
    return
end
error('Test_polyopvar_subsref failed (%s): expected an error containing "%s".',label,needle)
end
