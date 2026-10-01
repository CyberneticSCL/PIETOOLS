function test_index_maps()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_INDEX_MAPS checks the shared index maps of sopvar/misc/conventions,
%
%   gamma_of_cell, cell_of_gamma       parameter cell <-> multi-index gamma
%   multiindex_grid                    product grids, both orders
%   kron_strides, kron_split           Kronecker monomial positions
%   monomial_position_map/_position    degree -> position in one basis
%
% against two references:
%
% (A) the DEFINITIONS, by brute force and without any of the routines
%     under test: a cell index is the linear index MATLAB assigns to a
%     subscript of a 3 x ... x 3 array (read off an array holding its own
%     linear indices); a Kronecker position is where kron() of unit
%     vectors puts its 1; a product grid is ndgrid's; a degree position is
%     find(Z==d)-1;
%
% (B) every local form these replace, copied VERBATIM below from the files
%     at b3808f0f (function names prefixed by their file), on all cells for
%     0..5 shared variables and on random inputs in 1..4 directions: same
%     values, same orientation, same class. Where a local form errors or is
%     wrong on a degenerate input, the difference is asserted explicitly and
%     named in the helper's header.
%
% The copies are a frozen record of the pre-merge code: once the call sites
% adopt the helpers, (B) still states what they computed before.
%
% Initial coding MMP, 09/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(20260930);
nchk = 0;

% ============================================ gamma <-> linear cell index
for n3 = 0:5
    nc = 3^n3;
    % (A) The array holding its own linear indices, read by subscript.
    X = reshape(1:nc,[3*ones(1,n3),1,1]);
    G = gamma_of_cell((1:nc)',n3);
    assert(isequal(size(G),[nc,n3]) && isa(G,'double'),'gamma table size, n3 = %d',n3);
    g = ones(1,n3);                             % odometer over {1,2,3}^n3
    seen = false(nc,1);
    for r = 1:nc
        s = num2cell(g);
        v = X(s{:});                            % the definition of 'linear index'
        assert(isequal(G(v,:),g),'gamma_of_cell table row %d, n3 = %d',v,n3);
        assert(isequal(gamma_of_cell(v,n3),g),'gamma_of_cell(%d,%d)',v,n3);
        assert(isequal(cell_of_gamma(g),v),'cell_of_gamma, n3 = %d',n3);
        seen(v) = true;
        for t = 1:n3                            % advance the odometer
            if g(t)<3,  g(t) = g(t)+1;  break,  end
            g(t) = 1;
        end
        nchk = nchk+3;
    end
    assert(all(seen),'odometer missed a cell, n3 = %d',n3);
    kk = cell_of_gamma(G);
    assert(isequal(kk,(1:nc)'),'cell_of_gamma of the table, n3 = %d',n3);
    assert(isequal(gamma_of_cell(1:nc,n3),G),'row k input, n3 = %d',n3);   % row in, rows out
    nchk = nchk+2;

    % (B) Every local form, on every cell.
    for k = 1:nc
        gk = G(k,:);
        if n3>0
            % canonicalize_multiplier (2 sites), is_canonical_multiplier:
            % neither reaches this with n3 = 0 (both return first).
            assert(isequal(old_ind2sub_idiom(k,n3),gk),'ind2sub idiom k = %d n3 = %d',k,n3);
            nchk = nchk+1;
        end
        assert(isequal(old_adjmap_gamma(k,n3),gk),'canonical_adjoint_map/lpi_eq_sdopvar k = %d',k);
        assert(isequal(lpivar_cdopvar_gamma_of(k,n3),gk),'lpivar_cdopvar gamma_of k = %d',k);
        assert(isequal(eq_opts_sopvar_gamma_of(k,n3),gk),'eq_opts_sopvar gamma_of k = %d',k);
        assert(isequal(eq_opts_sopvar_lin_of(gk,n3),k),'eq_opts_sopvar lin_of k = %d',k);
        assert(isequal(degbalance_core_lin(gk,n3),k),'degbalance_core reach loop k = %d',k);
        assert(isequal(pi_gamma_index(k,n3),gk),'pi_gamma_index k = %d',k);
        if n3>0
            assert(isequal(old_sub2ind(gk,n3),k),'lpi_eq_sdopvar sub2ind k = %d',k);
            % the adjoint pairing of lpi_eq_sdopvar: 2 <-> 3
            adj = gk;   adj(gk==2) = 3;     adj(gk==3) = 2;
            assert(isequal(cell_of_gamma(adj),old_sub2ind(adj,n3)),'adjoint cell k = %d',k);
            nchk = nchk+2;
        end
        nchk = nchk+6;
    end
end
% Random multi-indices, 1..4 directions, many rows at once.
for trial = 1:200
    n3 = randi(4);      r = randi(20);
    Gr = randi(3,r,n3);
    kr = cell_of_gamma(Gr);
    assert(iscolumn(kr) && numel(kr)==r,'cell_of_gamma orientation');
    for i = 1:r
        c = num2cell(Gr(i,:));
        assert(kr(i)==sub2ind([3*ones(1,n3),1],c{:}),'random cell_of_gamma vs sub2ind');
    end
    assert(isequal(gamma_of_cell(kr,n3),Gr),'random gamma_of_cell');
    nchk = nchk+r+2;
end
% Refused inputs. The local forms do not refuse an out-of-range cell: the
% ind2sub idiom returns a last entry above 3 and the mod form wraps onto
% another cell. sub2ind refuses an out-of-range subscript, as does
% cell_of_gamma; the lin_of loops do not.
expect_error(@() gamma_of_cell(0,2),'integers in 1..3^n3');
expect_error(@() gamma_of_cell(10,2),'integers in 1..3^n3');
expect_error(@() gamma_of_cell(1.5,2),'integers in 1..3^n3');
expect_error(@() cell_of_gamma([0 1]),'should be 1, 2 or 3');
expect_error(@() cell_of_gamma([4 1]),'should be 1, 2 or 3');
expect_error(@() cell_of_gamma([1.5 1]),'should be 1, 2 or 3');
assert(isequal(old_ind2sub_idiom(10,2),[1 4]),'ind2sub idiom out of range');
assert(isequal(eq_opts_sopvar_gamma_of(10,2),[1 1]),'eq_opts gamma_of wraps');
expect_error(@() old_sub2ind([4 1],2),'');
assert(isequal(eq_opts_sopvar_lin_of([4 1],2),4),'lin_of does not refuse');
assert(isequal(cell_of_gamma(zeros(1,0)),1) && isequal(gamma_of_cell(1,0),zeros(1,0)),'n3 = 0');
assert(isequal(gamma_of_cell(zeros(0,1),2),zeros(0,2)),'empty k');
nchk = nchk+10;

% ================================================ multiindex_grid
for trial = 1:400
    N = randi(5)-1;                             % 0..4 directions
    vals = cell(1,N);
    for d = 1:N
        v = randi(9,1,randi(4))-1;              % 1..4 values, repeats allowed
        if rand<0.5,    v = v(:);   end         % any orientation
        vals{d} = v;
    end
    A  = multiindex_grid(vals);
    As = multiindex_grid(vals,'first_slowest');
    Af = multiindex_grid(vals,'first_fastest');
    % (A) ndgrid is the first-fastest product; the first-slowest one is the
    % first-fastest product over the reversed directions, reversed.
    assert(isequal(A,ndgrid_table(vals)) && isequal(Af,A),'first_fastest vs ndgrid, N = %d',N);
    assert(isequal(As,fliplr(ndgrid_table(fliplr(vals)))),'first_slowest vs ndgrid, N = %d',N);
    assert(isa(A,'double') && isa(As,'double'),'class');
    % (B) the local forms. Enumeration of alpha: row values as the callers
    % build them.
    rv = cellfun(@(x) reshape(x,1,[]),vals,'UniformOutput',false);
    assert(isequal(A,copquadvar_alpha_grid(vals)),'alpha_grid');
    assert(isequal(A,sopquadvar_alpha_all(rv)),'sopquadvar inline');
    if N>0
        assert(isequal(A,enum_alpha_kron(vals)),'enum_alpha / reach grid');
        nchk = nchk+1;
    end
    % Degree tables: the exponent vectors of a basis are columns.
    cv = cellfun(@(x) x(:),vals,'UniformOutput',false);
    assert(isequal(As,degbalance_core_degree_table(cv)),'degbalance_core degree_table');
    assert(isequal(As,monomial_gather_degree_table(cv)),'monomial_gather degree_table');
    assert(isequal(As,merge_monomial_product_full_degmat(cv)),'merge_monomial_product full_degmat');
    assert(isequal(As,kron_degree_oracle(cv)),'degree table vs kron of unit vectors');
    nchk = nchk+9;
end
% The eq_opts_sopvar / degbalance_core forms on their own inputs.
for n3 = 1:4
    for sv = 0:2^n3-1
        sepv = bitget(sv,1:n3)>0;
        vals = cell(1,n3);
        for k = 1:n3
            if sepv(k), vals{k} = [1,4]; else, vals{k} = [1,2,3]; end
        end
        assert(isequal(multiindex_grid(vals),eq_opts_sopvar_enum_alpha(n3,sepv)),'enum_alpha');
        assert(isequal(multiindex_grid(vals),degbalance_core_enum_alpha(n3,sepv)),'enum_alpha dbc');
        % reach: cells the diagonal of block a populates
        A = multiindex_grid(vals);
        for r = 1:size(A,1)
            a = A(r,:);
            sets = cell(1,n3);
            for k = 1:n3
                if a(k)==1, sets{k} = 1; else, sets{k} = [2 3]; end
            end
            C = cell_of_gamma(multiindex_grid(sets));
            assert(isequal(C,eq_opts_sopvar_reach(a,n3)),'reach eq_opts');
            assert(isequal(C,degbalance_core_reach(a,n3)),'reach degbalance_core');
            nchk = nchk+2;
        end
        nchk = nchk+2;
    end
end
% Degenerate inputs: the documented differences.
assert(isequal(multiindex_grid({}),zeros(1,0)) && isequal(multiindex_grid({},'first_slowest'),zeros(1,0)),'N = 0');
assert(isequal(copquadvar_alpha_grid({}),zeros(1,0)),'alpha_grid N = 0');
expect_error(@() eq_opts_sopvar_enum_alpha(0,false(1,0)),'');      % indexes vals{1}
assert(isequal(multiindex_grid({[1 2],[],3}),zeros(0,3)),'empty set');
expect_error(@() copquadvar_alpha_grid({[1 2],[],3}),'');          % 0/0 replication
assert(isequal(multiindex_grid({zeros(0,1),[1;2]},'first_slowest'),zeros(0,2)),'empty basis');
assert(isequal(monomial_gather_degree_table({zeros(0,1),[1;2]}),zeros(0,2)),'mg empty basis');
% degbalance_core's max(size(D,1),1) turns the 0 rows after an empty basis
% into 1, and MATLAB's concatenation then drops the empty column: a 2 x 1
% table for a basis with no monomial.
assert(isequal(degbalance_core_degree_table({zeros(0,1),[1;2]}),[1;2]),'dbc empty basis');
Zrow = {[0 1 2],[0;1]};                                             % a ROW basis
assert(isequal(multiindex_grid(Zrow,'first_slowest'),degbalance_core_degree_table(Zrow)),'row basis');
assert(~isequal(size(monomial_gather_degree_table(Zrow)),[6,2]),'monomial_gather widens a row basis');
expect_error(@() multiindex_grid({1},'kron'),'first_fastest');
nchk = nchk+11;

% ============================================ kron_strides / kron_split
for trial = 1:300
    N = randi(5)-1;
    nvec = randi(4,1,N);
    if rand<0.5,    nvec = nvec(:);     end
    s = kron_strides(nvec);
    assert(isequal(size(s),[1,N]) && isa(s,'double'),'kron_strides shape');
    assert(isequal(s,canonicalize_multiplier_strides_of(nvec)),'strides_of (cm)');
    assert(isequal(s,canonical_adjoint_map_strides_of(nvec)),'strides_of (cam)');
    assert(isequal(s,lpivar_cdopvar_strides(nvec)),'strides (lpivar)');
    NT = prod([nvec(:).',1]);
    t = (0:NT-1)';
    if rand<0.5,    t = t.';    end             % any shape in, a column table out
    a = kron_split(t,s);
    assert(isequal(size(a),[NT,N]),'kron_split shape');
    assert(isequal(a,canonicalize_multiplier_split_index(t,s)),'split_index (cm)');
    assert(isequal(a,canonical_adjoint_map_split_index(t,s)),'split_index (cam)');
    assert(isequal(a,copquadvar_onesided_split(t(:)+1,nvec(:).')),'int_onesided split');
    % (A) kron of unit vectors puts its 1 at the composite position.
    for i = 1:NT
        e = 1;
        for d = 1:N
            u = zeros(nvec(d),1);   u(a(i,d)+1) = 1;
            e = kron(e,u);
        end
        assert(find(e)==i,'kron position, row %d',i);
        assert(a(i,:)*s(:)==t(i),'inverse');
    end
    % Regrouping a subset of the directions (int_onesided's iL/iR/i3).
    dL = find(rand(1,N)<0.5);
    assert(isequal(a(:,dL)*kron_strides(nvec(dL)).',copquadvar_onesided_regroup(a,dL,nvec(:).')),'regroup');
    % tensor_index = multiindex_grid(idx)*st(:), idx subsets as rows
    idx = cell(1,N);
    for d = 1:N
        idx{d} = find(rand(1,nvec(d))<0.7)-1;
    end
    lin = multiindex_grid(idx)*s(:);
    assert(isequal(lin,lpivar_cdopvar_tensor_index(idx,s)),'tensor_index');
    nchk = nchk+NT*2+10;
end
assert(isequal(kron_strides([]),ones(1,0)) && isequal(kron_split(3,ones(1,0)),zeros(1,0)),'N = 0');
assert(isequal(multiindex_grid({})*ones(1,0).',lpivar_cdopvar_tensor_index({},ones(1,0))),'tensor_index N = 0');
nchk = nchk+2;

% ============================ monomial_position_map / monomial_position
for trial = 1:300
    n = randi(6);
    Z = randperm(10,n)'-1;                      % distinct exponents
    if rand<0.3,    Z = Z.';    end
    map = monomial_position_map(Z);
    assert(isequal(map,canonicalize_multiplier_degree_lookup(Z)),'degree_lookup (cm)');
    assert(isequal(map,canonical_adjoint_map_degree_lookup(Z)),'degree_lookup (cam)');
    deg = Z(randi(n,1,randi(8)));
    if rand<0.5,    deg = deg(:);   end
    out = monomial_position(map,deg);
    % (A) the position by search
    want = zeros(numel(deg),1);
    for i = 1:numel(deg),   want(i) = find(Z==deg(i))-1;   end
    assert(isequal(out,want),'monomial_position vs find');
    assert(isequal(out,canonicalize_multiplier_lookup(map,deg)),'lookup (cm)');
    assert(isequal(out,canonical_adjoint_map_lookup(map,deg)),'lookup (cam)');
    % An absent degree: same message as each copy.
    miss = setdiff(0:11,Z);     dm = [deg(:); miss(1)];
    m1 = errmsg(@() monomial_position(map,dm));
    m2 = errmsg(@() canonicalize_multiplier_lookup(map,dm));
    m3 = errmsg(@() monomial_position(map,dm,"Adjoint monomial basis"));
    m4 = errmsg(@() canonical_adjoint_map_lookup(map,dm));
    assert(~isempty(m1) && strcmp(m1,m2) && ~isempty(m3) && strcmp(m3,m4),'lookup messages');
    nchk = nchk+7;
end
assert(isequal(monomial_position_map([]),zeros(0,1)),'empty basis');
% Invalid exponents: the stricter copy's message; the other copy fails
% inside MATLAB on every one of them.
for bad = {[-1;0],[0;1.5],[2;-2],NaN}
    expect_error(@() monomial_position_map(bad{1}),'nonnegative integers');
    expect_error(@() canonicalize_multiplier_degree_lookup(bad{1}),'');
    nchk = nchk+2;
end

fprintf('test_index_maps: %d checks passed.\n',nchk);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Oracles from the definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = ndgrid_table(vals)
% Product grid, first direction fastest, by ndgrid.
N = numel(vals);
if N==0,    A = zeros(1,0);     return,     end
c = cell(1,N);
v = cellfun(@(x) reshape(x,[],1),vals,'UniformOutput',false);
if N==1
    c{1} = v{1};
else
    [c{:}] = ndgrid(v{:});
end
A = zeros(numel(c{1}),N);
for d = 1:N,    A(:,d) = c{d}(:);   end
end

function D = kron_degree_oracle(Z)
% Row t: the exponents of the t-th entry of kron(Z_1,...,Z_N), with the
% position of each factor found from kron() of unit vectors.
N = numel(Z);   n = cellfun(@numel,Z);  NT = prod([n,1]);
D = zeros(NT,N);
for d = 1:N
    % the entry t of kron(1,..,Z_d-index vector,..,1) is the Z_d position
    v = 1;
    for e = 1:N
        if e==d,    v = kron(v,(1:n(e))');
        else,       v = kron(v,ones(n(e),1));
        end
    end
    D(:,d) = reshape(Z{d}(v),[],1);
end
if NT==0,   D = zeros(0,N);     end
end

function expect_error(f,frag)
try
    f();
catch ME
    if ~isempty(frag) && ~contains(ME.message,frag)
        error('test_index_maps: expected ''%s'', got: %s',frag,ME.message);
    end
    return
end
error('test_index_maps: expected an error, none was raised.');
end

function m = errmsg(f)
m = '';
try
    f();
catch ME
    m = ME.message;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Local forms, VERBATIM at b3808f0f (renamed <file>_<function>)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gam = old_ind2sub_idiom(k,n3)
% canonicalize_multiplier (pass 1 and pass 3), is_canonical_multiplier
sz_C = [3*ones(1,n3),1];
    idcs = cell(1,n3);
    [idcs{:}] = ind2sub(sz_C,k);
    gam = cell2mat(idcs);
end

function gam = old_adjmap_gamma(k,n3)
% canonical_adjoint_map, lpi_eq_sdopvar
sz_C = [3*ones(1,n3),1];
    gam = ones(1,n3);
    if n3>0
        idcs = cell(1,n3);
        [idcs{:}] = ind2sub(sz_C,k);
        gam = cell2mat(idcs);
    end
end

function k = old_sub2ind(adj,n3)
% lpi_eq_sdopvar (adjoint cell of a self-adjoint program)
sz_C = [3*ones(1,n3),1];
        adj_cell = num2cell(adj);
        k = sub2ind(sz_C,adj_cell{:});
end

function k = degbalance_core_lin(g,n3)
% degbalance_core>reach, the inline loop, on row r = 1 of G = g
G = g;  r = 1;
    k = 1;
    for t = 1:n3, k = k + (G(r,t)-1)*3^(t-1); end
end

function gam = lpivar_cdopvar_gamma_of(g,n3)
% The gamma multi-index of linear cell g, first direction fastest, as the
% class's ind2sub over [3 ... 3] gives it.
if n3==0,   gam = zeros(1,0);   return,     end
c = cell(1,n3);
[c{:}] = ind2sub([3*ones(1,n3),1],g);
gam = cell2mat(c(1:n3));
end

function st = lpivar_cdopvar_strides(nvec)
% Stride of each variable in kron(Z{1},...,Z{N}), first variable slowest.
st = ones(1,numel(nvec));
for k = numel(nvec)-1:-1:1,     st(k) = st(k+1)*nvec(k+1);   end
end

function lin = lpivar_cdopvar_tensor_index(idx,st)
% Linear 0-based indices of the tensor product of per-variable index sets.
lin = 0;
for p = 1:numel(idx)
    lin = reshape(lin(:) + st(p)*reshape(idx{p},1,[]),[],1);
end
end

function g = eq_opts_sopvar_gamma_of(k,n3)
g = ones(1,n3);
for t = 1:n3
    g(t) = mod(floor((k-1)/3^(t-1)),3)+1;
end
end

function k = eq_opts_sopvar_lin_of(g,n3)
k = 1;
for t = 1:n3
    k = k + (g(t)-1)*3^(t-1);
end
end

function C = eq_opts_sopvar_reach(a,n3)
sets = cell(1,n3);
for k = 1:n3
    if a(k)==1
        sets{k} = 1;
    else
        sets{k} = [2 3];
    end
end
G = sets{1}(:);
for k = 2:n3
    G = [kron(ones(numel(sets{k}),1),G), ...
         kron(sets{k}(:),ones(size(G,1),1))];
end
C = zeros(size(G,1),1);
for r = 1:size(G,1)
    C(r) = eq_opts_sopvar_lin_of(G(r,:),n3);
end
end

function A = eq_opts_sopvar_enum_alpha(n3,sepv)
vals = cell(1,n3);
for k = 1:n3
    if sepv(k), vals{k} = [1,4]; else, vals{k} = [1,2,3]; end
end
A = vals{1}(:);
for k = 2:n3
    A = [kron(ones(numel(vals{k}),1),A), ...
         kron(vals{k}(:),ones(size(A,1),1))];
end
end

function A = enum_alpha_kron(vals)
% The kron construction of enum_alpha / reach, on given value sets.
n3 = numel(vals);
A = vals{1}(:);
for k = 2:n3
    A = [kron(ones(numel(vals{k}),1),A), ...
         kron(vals{k}(:),ones(size(A,1),1))];
end
end

function D = degbalance_core_degree_table(Z)
D = zeros(1,0);
for k = 1:numel(Z)
    D = [kron(D,ones(numel(Z{k}),1)), ...
         kron(ones(max(size(D,1),1),1),Z{k}(:))];
end
end

function C = degbalance_core_reach(a,n3)
sets = cell(1,n3);
for k = 1:n3
    if a(k)==1, sets{k} = 1; else, sets{k} = [2 3]; end
end
G = sets{1}(:);
for k = 2:n3
    G = [kron(ones(numel(sets{k}),1),G), ...
         kron(sets{k}(:),ones(size(G,1),1))];
end
C = zeros(size(G,1),1);
for r = 1:size(G,1)
    k = 1;
    for t = 1:n3, k = k + (G(r,t)-1)*3^(t-1); end
    C(r) = k;
end
end

function A = degbalance_core_enum_alpha(n3,sepv)
vals = cell(1,n3);
for k = 1:n3
    if sepv(k), vals{k} = [1,4]; else, vals{k} = [1,2,3]; end
end
A = vals{1}(:);
for k = 2:n3
    A = [kron(ones(numel(vals{k}),1),A), ...
         kron(vals{k}(:),ones(size(A,1),1))];
end
end

function D = monomial_gather_degree_table(Z)
D = zeros(1,0);
for k = 1:numel(Z)
    D = [kron(D,ones(size(Z{k}))), kron(ones(size(D,1),1),Z{k})];           %#ok<AGROW>
end
end

function degmat = merge_monomial_product_full_degmat(Zc)
degmat = zeros(1,0);
for i=1:numel(Zc)
    Zi = reshape(Zc{i},[],1);
    degmat = [kron(degmat,ones(numel(Zi),1)), kron(ones(size(degmat,1),1),Zi)];
end
end

function A = copquadvar_alpha_grid(vals)
nd = numel(vals);
szv = cellfun(@numel,vals);
nall = prod([szv,1]);
A = zeros(nall,nd);
rep = 1;
for d = 1:nd
    col = reshape(repmat(reshape(vals{d},1,[]),rep,1),[],1);
    A(:,d) = repmat(col,nall/(rep*szv(d)),1);
    rep = rep*szv(d);
end
end

function alpha_all = sopquadvar_alpha_all(vals)
% sopquadvar, the inline enumeration after 'vals' is built
n3 = numel(vals);
szv = cellfun(@numel,vals);                                                 % MMP, 09/12/2026
nall = prod([szv,1]);                                                       % MMP, 09/12/2026
alpha_all = zeros(nall,n3);                                                 % MMP, 09/12/2026
rep = 1;                                                                    % MMP, 09/12/2026
for k = 1:n3                                                                % MMP, 09/12/2026
    col = reshape(repmat(vals{k},rep,1),[],1);                              % MMP, 09/12/2026
    alpha_all(:,k) = repmat(col,nall/(rep*szv(k)),1);                       % MMP, 09/12/2026
    rep = rep*szv(k);                                                       % MMP, 09/12/2026
end                                                                         % MMP, 09/12/2026
end

function sub = copquadvar_onesided_split(tr,nout)
% copquadvar>int_onesided, the split of the transfer's output index
nv = numel(nout);
ns = numel(tr);
sub = zeros(ns,nv);
rem = tr-1;
for d = nv:-1:1
    sub(:,d) = mod(rem,nout(d));
    rem = floor(rem/nout(d));
end
end

function iL = copquadvar_onesided_regroup(sub,dL,nout)
% copquadvar>int_onesided, one group index from the split
ns = size(sub,1);
iL = zeros(ns,1);   for t = 1:numel(dL),  iL = iL*nout(dL(t))+sub(:,dL(t));  end
end

function s = canonicalize_multiplier_strides_of(nvec)
% Stride of each variable in kron(Z{1},...,Z{N}), first variable slowest.

N = numel(nvec);
s = ones(1,N);
for k = N-1:-1:1
    s(k) = s(k+1)*nvec(k+1);
end

end

function a = canonicalize_multiplier_split_index(idx,stride)
% Split a zero-based monomial index into zero-based per-variable indices.

N = numel(stride);
a = zeros(numel(idx),N);
rem = idx(:);
for k = 1:N
    a(:,k) = floor(rem/stride(k));
    rem = rem - a(:,k)*stride(k);
end

end

function map = canonicalize_multiplier_degree_lookup(Zvec)
% Position of a degree within a basis, as a zero-based index.

d = double(Zvec(:));
if isempty(d)
    map = zeros(0,1);
    return
end
map = -ones(max(d)+1,1);
map(d+1) = 0:numel(d)-1;

end

function out = canonicalize_multiplier_lookup(map,deg)
% Position of each requested degree in a basis, as a zero-based index.

deg = double(deg(:));
out = -ones(numel(deg),1);
in_range = deg>=0 & deg==round(deg) & deg+1<=numel(map);
out(in_range) = map(deg(in_range)+1);
bad = find(out<0,1);
if ~isempty(bad)
    error("Monomial basis does not contain degree "+num2str(deg(bad))...
          +"; the bases were not enlarged correctly.")
end

end

function s = canonical_adjoint_map_strides_of(nvec)
% Stride of each variable in kron(Z{1},...,Z{N}), first variable slowest.

N = numel(nvec);
s = ones(1,N);
for k = N-1:-1:1
    s(k) = s(k+1)*nvec(k+1);
end

end

function a = canonical_adjoint_map_split_index(idx,stride)
% Split a zero-based monomial index into zero-based per-variable indices.

N = numel(stride);
a = zeros(numel(idx),N);
rem = idx(:);
for k = 1:N
    a(:,k) = floor(rem/stride(k));
    rem = rem - a(:,k)*stride(k);
end

end

function map = canonical_adjoint_map_degree_lookup(Zvec)
% Position of a degree within a basis, as a zero-based index, stored as a
% direct lookup table over the shifted degree.

d = double(Zvec(:));
if isempty(d)
    map = zeros(0,1);
    return
end
if any(d<0) || any(d~=round(d))
    error("Monomial exponents should be nonnegative integers.")
end
map = -ones(max(d)+1,1);
map(d+1) = 0:numel(d)-1;

end

function out = canonical_adjoint_map_lookup(map,deg)
% Position of each requested degree in a basis, as a zero-based index.

deg = double(deg(:));
out = -ones(numel(deg),1);
in_range = deg>=0 & deg==round(deg) & deg+1<=numel(map);
out(in_range) = map(deg(in_range)+1);
bad = find(out<0,1);
if ~isempty(bad)
    error("Adjoint monomial basis does not contain degree "...
          +num2str(deg(bad))+"; the bases were not enlarged correctly.")
end

end
