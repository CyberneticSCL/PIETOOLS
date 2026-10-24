% Test_polyopvar_subsref.m -- verify @polyopvar/subsref.m, in two parts:
% (1) '()': row indexing, v(2,:), and the error for a column-restricted index;
% (2) '{}': f{i,j} must return a standalone 'polyopvar' object representing the term
%     f.C.ops{i,j} paired with monomial f.Z(j,:), i.e. the (i,j) contribution to f(x) =
%     C*Z(x). Summing apply_poly over every populated column j, for a fixed row i,
%     must reconstruct apply_poly(f,.) restricted to row i exactly. Also f{j,:}, the
%     nonempty terms of component j with their monomials, for f with diagonal blocks,
%     [y1;y2], [T1*y1;y1.*y2;y2], vector-valued blocks, y3, [y3;y1], [Top*y1;y2], rows
%     with several terms, and 1 and 2 spatial variables; f{:,k} (zero rows for the rows
%     without a block in column k, their dimensions taken from f.C). The values of the
%     blocks, xb{j}, come from the operators and polynomials directly (apply_sopvar),
%     independent of apply_poly(f) and of the class under test.

% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: '{}' with ':' and numeric column vectors, diagonal and
%   vector-valued blocks, multiple rows, linear index over blocks; f{:,k} for rows
%   without a block; '()' part merged from Test_polyopvar_cat_subsref (vertcat is in
%   Test_polyopvar_vertcat);

% NOTES:
% - '()' with an explicit column index (v(2,1)) is not supported: @polyopvar/subsref.m
%   forwards [row,col] to the 'tensopmat' subsref, whose columns are one block per
%   monomial, not per matrix column. Row-only v(2,:) is not affected. The test checks
%   that the column index errors with the message 'number of monomials must match', so
%   it fails once column indices are supported.
% - f{i,j} must return a 'tensopmat' coefficient (not the bare tensopvar/intop that
%   C.ops{i,j} unwraps to), else out.dims.out / out.depmat.out fail.
% - A (row,col) pair whose own block is empty and which is the sole entry of its
%   extracted row and column has no dimensions to extract: an error is expected.
% - 'end' evaluates to 1 for a 'polyopvar' (MATLAB's size of the scalar object): f{j,end}
%   is f{j,1}, f(end,:) is f(1,:). Not supported, and not tested except for this value.

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

% '()' subsref, row-only
v2 = v(2,:);
Vv2 = apply_poly(v2,{xv{1},yv{1}});
errv2 = cleanpoly(Vv2-yv{1},1e-10);
if max(max(abs(errv2.C)))>tol*max(1,max(max(abs(yv{1}.C))))
    error("Test_polyopvar_subsref failed: row-subsref value mismatch.")
end

% '()' subsref, genuine COLUMN restriction: must error (see NOTES)
f = x + x.*x;
threw = false;
try
    f1 = f(1,1); %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,'number of monomials must match')
        error("Test_polyopvar_subsref failed: column-restricted '()' errored with unexpected message: %s",ME.message)
    end
end
if ~threw
    error("Test_polyopvar_subsref failed: column-restricted '()' should have errored (still-open gap, see header note).")
end
end

function check_curly_subsref()
tol = 1e-8;
ntrials = 5;

pvar s1 s2
x = polyopvar('x',[s1,s2],[0,1]);
y = polyopvar('y',[s1,s2],[0,1]);
f = x + x.*x;       % 2 distinct monomial terms: bare x, and x^2
v = [f;y];          % row 1: f's 2 terms; row 2: y's single (different) term

[nr,nc] = size(v.C);
if nc<2
    error("Test_polyopvar_subsref failed: test setup expected >=2 monomial columns.")
end

for trial=1:ntrials

xv = rand_poly([1,1],{'s1','s2'},2,5);
yv = rand_poly([1,1],{'s1','s2'},2,5);
Yfull = apply_poly(v,{xv,yv});

for i=1:nr
    Ysum = 0;
    for j=1:nc
        op_ij = v.C.ops{i,j};
        if isempty(op_ij)
            continue
        end
        vij = v{i,j};
        if ~isa(vij,'polyopvar')
            error("Trial %d: v{%d,%d} did not return a 'polyopvar' object.",trial,i,j)
        elseif ~isequal(vij.Z,v.Z(j,:))
            error("Trial %d: v{%d,%d}.Z does not match v.Z(%d,:).",trial,i,j,j)
        end
        Ysum = Ysum + apply_poly(vij,{xv,yv});
    end
    erri = cleanpoly(Ysum-Yfull(i),1e-10);
    if max(max(abs(erri.C)))>tol*max(1,max(max(abs(Yfull(i).C))))
        error("Trial %d: sum of extracted terms for row %d does not reconstruct the original row.",trial,i)
    end
end

end

%%%% Linear indexing must match the corresponding (row,col) pair
v11 = v{1,1};   vlin1 = v{1};
if ~isequal(v11.Z,vlin1.Z) || ~isequal(v11.C.ops,vlin1.C.ops)
    error("Test_polyopvar_subsref failed: linear index v{1} does not match v{1,1}.")
end

%%%% Out-of-range indices must error
threw = false;
try
    bad = v{1,nc+1}; %#ok<NASGU>
catch
    threw = true;
end
if ~threw
    error("Test_polyopvar_subsref failed: out-of-range column index should have errored.")
end
threw = false;
try
    bad = v{nr+1,1}; %#ok<NASGU>
catch
    threw = true;
end
if ~threw
    error("Test_polyopvar_subsref failed: out-of-range row index should have errored.")
end

%%%% A (row,col) pair whose sole (1x1) block is empty cannot be extracted
empty_pair_found = false;
for i=1:nr
    for j=1:nc
        if isempty(v.C.ops{i,j})
            empty_pair_found = true;
            threw = false;
            try
                bad = v{i,j}; %#ok<NASGU>
            catch
                threw = true;
            end
            if ~threw
                error("Test_polyopvar_subsref failed: extracting an all-empty (%d,%d) block should have errored.",i,j)
            end
        end
    end
end
if ~empty_pair_found
    error("Test_polyopvar_subsref failed: test setup expected at least one empty (row,col) pair.")
end

%%%% ':' and numeric column vectors, diagonal and vector-valued blocks
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

        % v{j,:}: the nonempty terms of component j, with their monomials
        for j=1:nb
            vj = v{j,:};
            lbl = sprintf('%s, v{%d,:}',label,j);
            if ~isa(vj,'polyopvar') || size(vj,1)~=bs(j)
                error("Test_polyopvar_subsref failed (%s): wrong class or size.",lbl)
            elseif ~isequal(size(vj.C.ops),[1,nnz(pop(j,:))]) || any(cellfun(@isempty,vj.C.ops(:)))
                error("Test_polyopvar_subsref failed (%s): not the nonempty terms (%s).",lbl,mat2str(size(vj.C.ops)))
            elseif ~isequal(full(vj.Z),full(v.Z(pop(j,:),:)))
                error("Test_polyopvar_subsref failed (%s): monomials do not match those of the nonempty terms.",lbl)
            end
            assert_poly_equal(apply_poly(vj,xc),xb{j},tol,lbl);
            % an explicit list of all columns must give the same as ':'
            vj2 = v{j,1:nbc};
            if ~same_terms(vj,vj2)
                error("Test_polyopvar_subsref failed (%s): v{j,1:end} differs from v{j,:}.",lbl)
            end
            % a reversed list of the populated columns: same values, reversed monomials
            cj = find(pop(j,:));
            vj3 = v{j,fliplr(cj)};
            if ~isequal(full(vj3.Z),flipud(full(vj.Z)))
                error("Test_polyopvar_subsref failed (%s): column order of the index not kept.",lbl)
            end
            assert_poly_equal(apply_poly(vj3,xc),xb{j},tol,[lbl ', reversed columns']);
            % the terms one at a time sum to the component
            Ysum = 0;
            for jc=cj
                Ysum = Ysum + apply_poly(v{j,jc},xc);
            end
            assert_poly_equal(Ysum,xb{j},tol,[lbl ', sum of v{j,k}']);
            ncases = ncases+1;
        end

        % stacking v{j,:} over all j, and v{:,:}, reproduce v
        Vs = v{1,:};
        for j=2:nb,     Vs = [Vs;v{j,:}];    end                %#ok<AGROW>
        assert_poly_equal(apply_poly(Vs,xc),xv,tol,[label ', [v{1,:};v{2,:};...]']);
        Va = v{:,:};
        assert_poly_equal(apply_poly(Va,xc),xv,tol,[label ', v{:,:}']);
        if ~isequal(full(Va.Z),full(v.Z))
            error("Test_polyopvar_subsref failed (%s): v{:,:} does not retain all monomials.",label)
        end
        % several rows at once, in the given order
        if nb>=2
            for rows = {[1,nb],[nb,1],2:nb}
                r = rows{1};
                vr = v{r,:};
                assert_poly_equal(apply_poly(vr,xc),vertcat(xb{r}),tol,sprintf('%s, v{%s,:}',label,mat2str(r)));
            end
        end

        % term k of all rows, also if some rows have no block in column k (zero rows, dimensions
        % kept): the terms over k sum to the function
        Ysum = 0;
        for k=1:nbc
            vk = v{:,k};
            if size(vk,1)~=sum(bs) || size(vk.Z,1)~=1 || ~isequal(full(vk.Z),full(v.Z(k,:)))
                error("Test_polyopvar_subsref failed (%s): v{:,%d} is not term %d of all rows.",label,k,k)
            end
            Ysum = Ysum + apply_poly(vk,xc);
            ncases = ncases+1;
        end
        assert_poly_equal(Ysum,xv,tol,[label ', sum of v{:,k}']);

        % linear index over the blocks of C.ops
        for k=1:nbr*nbc
            [r,c] = ind2sub([nbr,nbc],k);
            if pop(r,c)
                vk = v{k};      vrc = v{r,c};
                if ~same_terms(vk,vrc)
                    error("Test_polyopvar_subsref failed (%s): linear index v{%d} differs from v{%d,%d}.",label,k,r,c)
                end
            else
                assert_error(@() v{k},'no nonempty term',sprintf('%s, linear index to an empty block',label));
            end
        end
        % indices beyond the blocks (the scalar dimensions size(v.C) are larger for vector-valued blocks)
        assert_error(@() v{nbr+1,1},'Row indices cannot exceed',[label ', row beyond the blocks']);
        assert_error(@() v{1,nbc+1},'Column indices cannot exceed',[label ', column beyond the blocks']);
        ncases = ncases+3;
    end
end

%%%% Empty selections, and the indices that remain unsupported
vv = [ya;yb];
assert_error(@() vv{':'},'Linear indexing',"linear ':'");
assert_error(@() vv{1,1,1},'At most two indices',"three indices");
assert_error(@() vv{[],:},'no nonempty term',"empty row selection");
assert_error(@() vv{1,[]},'no nonempty term',"empty column selection");
assert_error(@() vv{1,2},'no nonempty term',"an empty (row,col) pair, 2x2 blocks");
% subsequent references are applied to the extracted object
if ~isequal(full(vv{2,:}.Z),2) || ~isequal(vv{2,:}.dims.out,1)
    error("Test_polyopvar_subsref failed: chained reference after vv{2,:}.")
end
% 'end' evaluates to 1, see NOTES
if ~same_terms(vv{1,end},vv{1,1}) || ~same_terms(vv{end,1},vv{1,1})
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
